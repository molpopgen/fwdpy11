#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <fwdpy11/genetic_values/MultivariateGSSmo.hpp>

namespace py = pybind11;


void
init_MultivariateGSSmo(py::module& m)
{
    py::class_<fwdpy11::MultivariateGSSmo,
               fwdpy11::MultivariateGeneticValueToFitnessMap>(
        m, "MultivariateGSSmo",
        "Multivariate Gaussian stabilizing selection with moving optima.")
        .def(py::init([](py::array_t<std::uint32_t> timepoints,
                         py::array_t<double> optima, double VS) {
                 auto t = timepoints.unchecked<1>();
                 auto o = optima.unchecked<2>();

                 std::vector<std::uint32_t> it(t.data(0),
                                               t.data(0) + t.shape(0));
                 std::vector<double> io(
                     o.data(0, 0), o.data(0, 0) + o.shape(0) * o.shape(1));

                 return fwdpy11::MultivariateGSSmo(std::move(it),
                                                   std::move(io), VS);
             }),
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
        )delim")
        .def(py::pickle(
            [](const fwdpy11::MultivariateGSSmo& self) {
                return self.pickle();
            },
            [](py::object o) {
                auto t = o.cast<py::tuple>();
                auto l = t[0].cast<py::list>();
                std::vector<std::uint32_t> tp;
                for (auto i : l)
                    {
                        tp.push_back(i.cast<std::uint32_t>());
                    }
                l = t[1].cast<py::list>();
                std::vector<double> optima;
                for (auto i : l)
                    {
                        optima.push_back(i.cast<double>());
                    }
                double vs = t[2].cast<double>();
                return fwdpy11::MultivariateGSSmo(std::move(tp),
                                                  std::move(optima), vs);
            }))
        .def("__eq__", [](const fwdpy11::MultivariateGSSmo& lhs,
                          const fwdpy11::MultivariateGSSmo& rhs) {
            return lhs.timepoints == rhs.timepoints && lhs.optima == rhs.optima
                   && lhs.VS == rhs.VS;
        });
}
