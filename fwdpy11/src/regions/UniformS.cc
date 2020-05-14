#include <pybind11/pybind11.h>
#include <fwdpy11/regions/UniformS.hpp>

namespace py = pybind11;

void
init_UniformS(py::module& m)
{
    py::class_<fwdpy11::UniformS, fwdpy11::Sregion>(m, "_ll_UniformS")
        .def(py::init([](double beg, double end, double weight, double lo, double hi,
                         double h, bool coupled, std::uint16_t label, double scaling) {
                 return fwdpy11::UniformS(
                     fwdpy11::Region(beg, end, weight, coupled, label), scaling, lo, hi,
                     h);
             }),
             py::arg("beg"), py::arg("end"), py::arg("weight"), py::arg("lo"),
             py::arg("hi"), py::arg("h"), py::arg("coupled"), py::arg("label"),
             py::arg("scaling"));
}

