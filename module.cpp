#include "nbody.h"

PYBIND11_MODULE(nbody, m)
{
m.doc() = "nah";

py::class_<SimualationParameters>(m, "SimulationParameters")
        .def(py::init<>())
        .def_readwrite("G", &SimualationParameters::G)
        .def_readwrite("easing", &SimualationParameters::easing)
        .def_readwrite("friction", &SimualationParameters::friction)
        .def_readwrite("dt", &SimualationParameters::dt);

py::class_<EulerSimulation>(m, "EulerSimulation")
        .def(py::init<SimualationParameters&, Mref, Mref, Mref>())
        .def("step", &EulerSimulation::step)
        ;



//test if pybind works
m.def("incr_matrix_any", [](py::EigenDRef<Eigen::MatrixXd> m, double v) {
m += Eigen::MatrixXd::Constant(m.rows(), m.cols(), v);
return m;
}, py::return_value_policy::reference);

//helpers
m.def("print_matrix", &print_matrix);
m.def("init_threads", &threading);
}


//junkyard
//
//py::class_<Integrator>(m, "Integrator");
//py::class_<EulerIntegrator,Integrator>(m, "EulerIntegrator")
//.def(py::init<>());