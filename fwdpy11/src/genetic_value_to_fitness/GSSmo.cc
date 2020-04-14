#include <tuple>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <fwdpy11/genetic_value_to_fitness/GSSmo.hpp>

namespace py = pybind11;

static const auto INIT_TUPLES =
    R"delim(
:param optima: Model parameters over time
:type optima: list

Each element of optima must be a tuple of 
(generation, optimal trait value, VS)

.. deprecated:: 0.7.1
)delim";

static const auto INIT_OPTIMA =
    R"delim(
:param optima: Optima, sorted by time
:type optima: list

Each element of `optima` must be an instance
of :class:`fwdpy11.Optimum`.

.. versionadded:: 0.7.1
)delim";

void
init_GSSmo(py::module& m)
{
    py::class_<fwdpy11::GSSmo, fwdpy11::GeneticValueIsTrait>(
        m, "GSSmo", "Gaussian stabilizing selection with a moving optimum.")
        .def(py::init([](std::vector<std::tuple<std::uint32_t, double, double>> tuples) {
                 PyErr_WarnEx(PyExc_DeprecationWarning,
                              "Construction with tuples is deprecated.  Please use list "
                              "of fwdpy11.Optimum instead",
                              0);
                 std::vector<fwdpy11::Optimum> optima;
                 for (auto& t : tuples)
                     {
                         optima.emplace_back(std::get<0>(t), std::get<1>(t),
                                             std::get<2>(t));
                     }
                 return fwdpy11::GSSmo(std::move(optima));
             }),
             py::arg("optima"), INIT_TUPLES)
        .def(py::init<std::vector<fwdpy11::Optimum>>(), py::arg("optima"), INIT_OPTIMA)
        .def_readonly("VS", &fwdpy11::GSSmo::VS)
        .def_readonly("optimum", &fwdpy11::GSSmo::opt)
        .def_property_readonly(
            "opt",
            [](const fwdpy11::GSSmo& self) {
                PyErr_WarnEx(PyExc_DeprecationWarning,
                             "GSSmo.opt is deprecated.  Use GSSmo.optimum instead.", 0);
                return self.opt;
            })
        .def_readonly("optima", &fwdpy11::GSSmo::optima)
        .def(py::pickle([](const fwdpy11::GSSmo& g) { return g.pickle(); },
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
