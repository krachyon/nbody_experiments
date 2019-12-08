#include <vector>
#include <cmath>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <iostream>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-copy"

#include <omp.h>

#include <Eigen/LU>

#pragma GCC diagnostic pop

namespace py = pybind11;

using VectorType = Eigen::VectorXd;
using MatrixType = Eigen::MatrixXd;
using Mref = py::EigenDRef<MatrixType>;
using Vref = py::EigenDRef<VectorType>;


struct SimualationParameters{
    SimualationParameters() = default;
    double G = 0.1;
    double easing = 0.1;
    double dt = 0.001;
    double friction = 0;
};

struct Simulation
{
    virtual void step(size_t) = 0;
};

struct EulerSimulation: public Simulation
{
    EulerSimulation(SimualationParameters& params, Mref positions, Mref velocities, Mref masses):
            params(params), xs(positions), vs(velocities), masses(masses)
            , _dv{xs.rows(), 2}
    {
        assert(xs.rows() == vs.rows() && masses.rows() == vs.rows());
    }

    EulerSimulation() = delete;

    inline MatrixType& accel(Mref xs)
    {
        _dv.setZero();
        for(Eigen::Index i=0; i!= xs.rows(); ++i) {
            for(Eigen::Index j=i; j!= xs.rows(); ++j)
            {
                if(i==j)
                    continue;
                auto const diff = xs.row(j)-xs.row(i);
                auto const norm = diff.norm();
                auto da = params.G * masses[j] * diff/ (norm*norm*norm + params.easing);
                _dv.row(i) += da;
                _dv.row(j) -= da;
            }
        }
        return _dv;
    }


    void step(size_t iters) override
    {
        for (size_t it = 0; it!=iters; ++it) {
            vs += accel(xs) * params.dt;
            vs -= vs*params.friction;
            xs += params.dt * vs;
        }

    }
    SimualationParameters params;
    Mref xs;
    Mref vs;
    Vref masses;
    MatrixType _dv;

};

struct VerletSimulation: public Simulation
{
    //TODO broken
    // y'' = f(x,y)
    //y_i+1 = 2y_i − y_i−1 + ∆x^2 f (x, y).

    VerletSimulation(SimualationParameters params, Mref xs, Mref vs , Vref mass)
    : params(params), xs(xs), masses(mass), A{xs.rows(), 2}
    {
        assert(xs.rows() == vs.rows() && masses.rows() == vs.rows());

        //first of all, a single step with different integration to fill out xs_2
        xs_prev = xs;
        EulerSimulation(params, xs, vs, mass).step(1);
    }


    void step(size_t iters) override {
        A.setZero();
        for (size_t _ = 0; _ != iters; ++_) {
            for (Eigen::Index i = 0; i != xs.rows(); ++i) {
                for (Eigen::Index j = 0; j != xs.rows(); ++j) {
                    if (j==i)
                        continue;
                    // particle i feels F,  j feels -F
                    VectorType diff = xs.row(i) - xs.row(j);
                    diff = diff / pow(diff.norm(),3);
                    auto dA = -1*params.G * masses[j] * diff;
                    A.row(i) += dA;
                    //F.row(j) -= dF;
                }
            }
            //verlet formula
            MatrixType xs_next = 2 * xs_prev - xs + pow(params.dt,2) * A;

            xs_prev = xs;
            xs = xs_next;
        }
    }


    SimualationParameters params;
    Mref xs;
    Vref masses;
    MatrixType xs_prev;
    MatrixType A;
};

struct LeapfrogSimulation: public Simulation
{
    LeapfrogSimulation(SimualationParameters params, Mref xs, Mref vs , Vref mass)
    : params(params), xs(xs), vs(vs), masses(mass),
        vs_half{vs.rows(), vs.cols()}, dv{vs.rows(),vs.cols()}
    {
        assert(xs.rows() == vs.rows() && masses.rows() == vs.rows());

        vs_half = vs + accel(xs) * params.dt / 2 ;
    }

    inline MatrixType& accel(Mref xs)
    {
        dv.setZero();
        for(Eigen::Index i=0; i!= xs.rows(); ++i) {
            for(Eigen::Index j=i; j!= xs.rows(); ++j)
            {
                if(i==j)
                    continue;
                auto const diff = xs.row(j)-xs.row(i);
                auto const norm = diff.norm();
                auto da = params.G * masses[j] * diff/ (norm*norm*norm + params.easing);
                dv.row(i) += da;
                dv.row(j) -= da;
            }
        }
        return dv;
    }

    void step(size_t iters) override
    {
        auto const dt = params.dt;

        for(size_t _=0;_!=iters;++_) {
//            vs_half = vs + accel(xs) * dt / 2 ;
//            xs += vs_half * params.dt;
//            vs = vs_half + accel(xs) * dt / 2;
//
//            vs -= vs_half*params.friction;
            xs += vs_half * dt;
            vs_half += accel(xs)*dt;
            vs_half -= vs_half * params.friction;
        }
    }

    SimualationParameters params;

    Mref xs;
    Mref vs;
    Vref masses;
    MatrixType vs_half;
    MatrixType dv;
};


void threading(size_t nthreads = 10)
{
    omp_set_num_threads(nthreads);
    Eigen::setNbThreads(nthreads);
    Eigen::initParallel();
}

void print_matrix(py::EigenDRef<MatrixType> mat)
{
    std::cout << mat << std::endl;
}
