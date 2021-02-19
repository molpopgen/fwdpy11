#include <pybind11/pybind11.h>
#include <fwdpy11/regions/GaussianS.hpp>
#include <fwdpy11/mutation_dominance/MutationDominance.hpp>

namespace py = pybind11;

void
init_GaussianS(py::module& m)
{
    py::class_<fwdpy11::GaussianS, fwdpy11::Sregion>(m, "_ll_GaussianS")
        .def(py::init([](double beg, double end, double weight, double sd, double h,
                         bool coupled, std::uint16_t label, double scaling) {
                 return fwdpy11::GaussianS(
                     fwdpy11::Region(beg, end, weight, coupled, label), scaling, sd, h);
             }),
             py::arg("beg"), py::arg("end"), py::arg("weight"), py::arg("sd"),
             py::arg("h"), py::arg("coupled"), py::arg("label"), py::arg("scaling"))
        .def(py::init([](double beg, double end, double weight, double sd,
                         const fwdpy11::MutationDominance& h, bool coupled,
                         std::uint16_t label, double scaling) {
                 return fwdpy11::GaussianS(
                     fwdpy11::Region(beg, end, weight, coupled, label), scaling, sd, h);
             }),
             py::arg("beg"), py::arg("end"), py::arg("weight"), py::arg("sd"),
             py::arg("h"), py::arg("coupled"), py::arg("label"), py::arg("scaling"));
}

