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
        m, "MultivariateGSSmo")
        .def(py::init([](py::array_t<std::uint32_t> timepoints,
                         py::array_t<double> optima, double VS) {
            auto t = timepoints.unchecked<1>();
            auto o = optima.unchecked<2>();

            std::vector<std::uint32_t> it(t.data(0), t.data(0) + t.shape(0));
            std::vector<double> io(o.data(0, 0),
                                   o.data(0, 0) + o.shape(0) * o.shape(1));

            return fwdpy11::MultivariateGSSmo(std::move(it), std::move(io),
                                              VS);
        }))
        .def(py::pickle(
            [](const fwdpy11::MultivariateGSSmo& self) {
                return self.pickle();
            },
            [](py::object o) {
                auto t = o.cast<py::tuple>();
                auto tp = t[0].cast<std::vector<std::uint32_t>>();
                auto optima = t[1].cast<std::vector<double>>();
                double vs = t[2].cast<double>();
                return fwdpy11::MultivariateGSSmo(std::move(tp),
                                                  std::move(optima), vs);
            }));
}
