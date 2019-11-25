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

struct EulerSimulation
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


    virtual void step(size_t iters)
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
