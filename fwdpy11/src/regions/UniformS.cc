#include <pybind11/pybind11.h>
#include <fwdpy11/regions/UniformS.hpp>

namespace py = pybind11;

static const auto INIT_DOCSTRING =
    R"delim(
Constructor

:param beg: the beginning of the region
:type beg: float
:param end: the end of the region
:type end: float
:param weight: the weight to assign
:type weight: float
:param lo: lower bound on s
:type lo: float
:param hi: upper bound on s
:type hi: float
:param h: the dominance
:type h: float
:param coupled: if True, the weight is converted to (end-beg)*weight
:type coupled: bool
:param label: Fill :attr:`fwdpy11.Mutation.label` with this value.
:type label: numpy.uint16
:param scaling: The scaling of the DFE
:type scaling: float

When coupled is True, the "weight" may be
interpreted as a "per base pair"
(or per unit, generally speaking) term.

Example:

.. testcode::

    #A simple case
    import fwdpy11
    #s is uniform on [-1, 0]
    uniformS = fwdpy11.UniformS(0,1,1,-1,0,0)
)delim";

void
init_UniformS(py::module& m)
{
    py::class_<fwdpy11::UniformS, fwdpy11::Sregion>(
        m, "UniformS", "Uniform distrubution of effect sizes")
        .def(py::init([](double beg, double end, double weight, double lo, double hi,
                         double h, bool coupled, std::uint16_t label, double scaling) {
                 return fwdpy11::UniformS(
                     fwdpy11::Region(beg, end, weight, coupled, label), scaling, lo, hi,
                     h);
             }),
             py::arg("beg"), py::arg("end"), py::arg("weight"), py::arg("lo"),
             py::arg("hi"), py::arg("h") = 1.0, py::arg("coupled") = true,
             py::arg("label") = 0, py::arg("scaling") = 1.0, INIT_DOCSTRING)
        .def_readonly("lo", &fwdpy11::UniformS::lo)
        .def_readonly("hi", &fwdpy11::UniformS::hi)
        .def_readonly("h", &fwdpy11::UniformS::dominance)
        .def("__repr__", &fwdpy11::UniformS::repr)
        .def(py::pickle([](const fwdpy11::UniformS& self) { return self.pickle(); },
                        [](py::tuple t) { return fwdpy11::UniformS::unpickle(t); }));
}

