#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <fwdpy11/genetic_value_to_fitness/MultivariateGSSmo.hpp>

namespace py = pybind11;

static const auto INIT_OPTIMA =
    R"delim(
:param optima: List of :class:`fwdpy11.PleiotropicOptima`
:type optima: list

.. versionadded:: 0.7.1
)delim";

void
init_MultivariateGSSmo(py::module& m)
{
    py::class_<fwdpy11::MultivariateGSSmo, fwdpy11::GeneticValueIsTrait>(
        m, "MultivariateGSSmo",
        "Multivariate Gaussian stabilizing selection with moving optima.")
        .def(py::init<const std::vector<fwdpy11::PleiotropicOptima>&>(), INIT_OPTIMA)
        .def(py::pickle(
            [](const fwdpy11::MultivariateGSSmo& self) { return self.pickle(); },
            [](py::object o) {
                auto l = o.cast<py::list>();
                std::vector<fwdpy11::PleiotropicOptima> optima;
                for (auto i : l)
                    {
                        auto j = i.cast<py::tuple>();
                        optima.emplace_back(j[0].cast<std::uint32_t>(),
                                            j[1].cast<std::vector<double>>(),
                                            j[2].cast<double>());
                    }
                return fwdpy11::MultivariateGSSmo(optima);
            }))
        .def("__eq__", [](const fwdpy11::MultivariateGSSmo& lhs,
                          const fwdpy11::MultivariateGSSmo& rhs) {
            return lhs.optima == rhs.optima;
        });
}
