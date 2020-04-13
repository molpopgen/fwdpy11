#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <fwdpy11/genetic_value_to_fitness/MultivariateGSSmo.hpp>

namespace py = pybind11;

static const auto INIT_DEPRECATED =
    R"delim(
:param timepoints: Time when the optima change
:type timepoints: numpy.array
:param optima: The optima corresponding to each time point
:type optima: numpy.ndarray
:param VS: Strength of stabilizing selection
:type VS: float

The rows of optima should correpond to the optimal trait values at
each time point.  The first row must correspond to a time point of zero.

The following example changes the optimum from 0 to 1 for the first
trait at 10N generations into a simulation:

.. testcode::

import fwdpy11
import numpy as np
popsize = 1000
timepoints = np.array([0,10*popsize], dtype=np.uint32)
ntraits = 3
optima = np.array(np.zeros(2*ntraits)).reshape((len(timepoints), ntraits))
optima[1,0] = 1
mvgssmo = fwdpy11.MultivariateGSSmo(timepoints,optima, 1)

.. deprecated:: 0.7.1
)delim";

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
        .def(py::init([](py::array_t<std::uint32_t> timepoints,
                         py::array_t<double> optima, double VS) {
                 PyErr_WarnEx(PyExc_DeprecationWarning,
                              "This __init__ function is deprecated.  Please use list "
                              "of fwdpy11.PleiotropicOptima instead",
                              0);

                 auto t = timepoints.unchecked<1>();
                 auto o = optima.unchecked<2>();

                 std::vector<fwdpy11::PleiotropicOptima> po;
                 for (py::ssize_t i = 0; i < t.size(); ++i)
                     {
                         std::vector<double> temp;
                         for (py::ssize_t j = 0; j < o.shape(1); ++j)
                             {
                                 temp.push_back(o(i, j));
                             }
                         po.emplace_back(t(i), std::move(temp), VS);
                     }
                 return fwdpy11::MultivariateGSSmo(po);
             }),
             INIT_DEPRECATED)
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
