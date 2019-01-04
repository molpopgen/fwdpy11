#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <fwdpy11/genetic_values/MultivariateGSS.hpp>

namespace py = pybind11;

void
init_MultivariateGSS(py::module& m)
{
    py::class_<fwdpy11::MultivariateGSS,
               fwdpy11::MultivariateGeneticValueToFitnessMap>(
        m, "MultivariateGSS")
        .def(py::init([](py::array_t<double> optima, double VS) {
            auto r = optima.unchecked<1>();
            std::vector<double> voptima(r.data(0), r.data(0) + r.shape(0));
            return fwdpy11::MultivariateGSS(std::move(voptima), VS);
        }));
}
