#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <fwdpy11/genetic_value_to_fitness/MultivariateGSS.hpp>

namespace py = pybind11;

void
init_MultivariateGSS(py::module& m)
{
    py::class_<fwdpy11::MultivariateGSS, fwdpy11::GeneticValueIsTrait>(
        m, "_ll_MultivariateGSS")
        .def(py::init([](py::array_t<double> optima, double VS) {
                 auto r = optima.unchecked<1>();
                 std::vector<double> voptima(r.data(0), r.data(0) + r.shape(0));
                 return fwdpy11::MultivariateGSS(
                     fwdpy11::PleiotropicOptima(std::move(voptima), VS));
             }),
             py::arg("optima"), py::arg("VS"))
        .def(py::init<const fwdpy11::PleiotropicOptima&>(), py::arg("optima"));
}
