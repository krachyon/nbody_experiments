#include <vector>

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
using Mref = py::EigenDRef<Eigen::MatrixXd>;
using Vref = py::EigenDRef<Eigen::VectorXd>;


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
            , _dv{xs.rows(), 2}, _dx{xs.rows(), 2}
    {
        assert(xs.rows() == vs.rows() && masses.rows() == vs.rows());
    }

    EulerSimulation() = delete;

    inline Eigen::Vector2d accel(size_t idx) const
    {
        Eigen::Vector2d dv{0., 0.};
        Eigen::Vector2d x = xs.row(idx);
        Eigen::Vector2d diff{0., 0.};


        size_t rows = xs.rows();
        //omp massively decreases performance here. WTF?
        //#pragma omp parallel for
        for (Eigen::Index i = 0; i!=rows; ++i) {
            if (i==idx)
                continue;

            Eigen::Vector2d const& otherx = xs.row(i);
            diff = otherx-x;
            dv += diff / diff.norm() * params.G * masses[i] / (diff.dot(diff) + params.easing);
        }

        return dv;
    }


    void step(size_t iters) override
    {
        for (size_t it = 0; it!=iters; ++it) {
            for (Eigen::Index i = 0; i != xs.rows(); ++i) {
                _dv.row(i) = accel(i);
                _dx.row(i) = (params.dt * vs.row(i));
            }
            vs += _dv;
            vs -= vs*params.friction;
            xs += _dx;
        }

    }
    SimualationParameters params;
    Mref xs;
    Mref vs;
    Vref masses;
    Eigen::MatrixXd _dv;
    Eigen::MatrixXd _dx;

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
                    Eigen::VectorXd diff = xs.row(i) - xs.row(j);
                    diff = diff / pow(diff.norm(),3);
                    auto dA = -1*params.G * masses[j] * diff;
                    A.row(i) += dA;
                    //F.row(j) -= dF;
                }
            }
            //verlet formula
            Eigen::MatrixXd xs_next = 2 * xs_prev - xs + pow(params.dt,2) * A;

            xs_prev = xs;
            xs = xs_next;
        }
    }


    SimualationParameters params;
    Mref xs;
    Vref masses;
    Eigen::MatrixXd xs_prev;
    Eigen::MatrixXd A;
};


void threading(size_t nthreads = 10)
{
    omp_set_num_threads(nthreads);
    Eigen::setNbThreads(nthreads);
    Eigen::initParallel();
}

void print_matrix(py::EigenDRef<Eigen::MatrixXd> mat)
{
    std::cout << mat << std::endl;
}
