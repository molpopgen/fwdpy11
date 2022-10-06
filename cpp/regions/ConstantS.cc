#include <pybind11/pybind11.h>
#include <fwdpy11/regions/ConstantS.hpp>
#include <fwdpy11/mutation_dominance/MutationDominance.hpp>

namespace py = pybind11;

void
init_ConstantS(py::module& m)
{
    py::class_<fwdpy11::ConstantS, fwdpy11::Sregion>(m, "_ll_ConstantS")
        .def(py::init([](double beg, double end, double weight, double s, double h,
                         bool coupled, std::uint16_t label, double scaling) {
                 return fwdpy11::ConstantS(
                     fwdpy11::Region(beg, end, weight, coupled, label), scaling, s, h);
             }),
             py::arg("beg"), py::arg("end"), py::arg("weight"), py::arg("s"),
             py::arg("h"), py::arg("coupled"), py::arg("label"), py::arg("scaling"))
        .def_readonly("esize", &fwdpy11::ConstantS::esize)
        .def(py::init([](double beg, double end, double weight, double s,
                         const fwdpy11::MutationDominance& h, bool coupled,
                         std::uint16_t label, double scaling) {
                 return fwdpy11::ConstantS(
                     fwdpy11::Region(beg, end, weight, coupled, label), scaling, s, h);
             }),
             py::arg("beg"), py::arg("end"), py::arg("weight"), py::arg("s"),
             py::arg("h"), py::arg("coupled"), py::arg("label"), py::arg("scaling"))
        .def_readonly("esize", &fwdpy11::ConstantS::esize);
}

