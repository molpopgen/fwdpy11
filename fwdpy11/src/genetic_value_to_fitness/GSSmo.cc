#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <fwdpy11/genetic_values/GeneticValueToFitness.hpp>

namespace py = pybind11;

void
init_GSSmo(py::module& m)
{
    py::class_<fwdpy11::GSSmo, fwdpy11::GeneticValueIsTrait>(
        m, "GSSmo", "Gaussian stabilizing selection with a moving optimum.")
        .def(
            py::init<std::vector<std::tuple<std::uint32_t, double, double>>>(),
            py::arg("optima"),
            R"delim(
            :param optima: Model parameters over time
            :type optima: list
            
            Each element of optima must be a tuple of 
            (generation, optimal trait value, VS)
            )delim")
        .def_readonly("VS", &fwdpy11::GSSmo::VS)
        .def_readonly("opt", &fwdpy11::GSSmo::opt)
        .def_readonly("optima", &fwdpy11::GSSmo::optima)
        .def(py::pickle(
            [](const fwdpy11::GSSmo& g) { return g.pickle(); },
            [](py::object o) {
                py::tuple t(o);
                if (t.size() != 4)
                    {
                        throw std::runtime_error("invalid object state");
                    }
                auto opt = t[0].cast<double>();
                auto VS = t[1].cast<double>();
                auto co
                    = t[2].cast<decltype(fwdpy11::GSSmo::current_optimum)>();
                auto optima = t[3].cast<decltype(fwdpy11::GSSmo::optima)>();
                auto rv = fwdpy11::GSSmo(std::move(optima));
                rv.opt = opt;
                rv.VS = VS;
                rv.current_optimum = co;
                return rv;
            }));
}
