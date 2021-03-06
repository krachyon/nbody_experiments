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
        .def_readwrite("xs", &EulerSimulation::xs)
        .def_readwrite("ys", &EulerSimulation::vs)
        .def_readwrite("ms", &EulerSimulation::masses)
        .def_readwrite("param", &EulerSimulation::params)
        ;

py::class_<VerletSimulation>(m, "VerletSimulation")
        .def(py::init<SimualationParameters&, Mref, Mref, Vref>())
        .def("step", &VerletSimulation::step)
        .def_readwrite("param", &VerletSimulation::params)
        .def_readwrite("xs", &VerletSimulation::xs)
        .def_readwrite("xs_prev", &VerletSimulation::xs_prev)
        .def_readwrite("A", &VerletSimulation::A)
        ;

py::class_<LeapfrogSimulation>(m, "LeapfrogSimulation")
        .def(py::init<SimualationParameters&, Mref, Mref, Vref>())
        .def("step", &LeapfrogSimulation::step)
        .def_readwrite("param", &LeapfrogSimulation::params)
        .def_readwrite("xs", &LeapfrogSimulation::xs)
        .def_readwrite("vs", &LeapfrogSimulation::vs)
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