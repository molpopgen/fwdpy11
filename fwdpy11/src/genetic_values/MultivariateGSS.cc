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
        }))
        .def(py::pickle(
            [](const fwdpy11::MultivariateGSS& self) { return self.pickle(); },
            [](py::object o) {
                py::tuple t = o.cast<py::tuple>();
                if (t.size() != 2)
                    {
                        throw std::invalid_argument("incorrect tuple size");
                    }
                std::vector<double> optima = t[0].cast<std::vector<double>>();
                double VS = t[1].cast<double>();
                return fwdpy11::MultivariateGSS(std::move(optima), VS);
            }));
}
