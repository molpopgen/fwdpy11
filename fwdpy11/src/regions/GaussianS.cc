#include <pybind11/pybind11.h>
#include <fwdpy11/regions/GaussianS.hpp>

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
:param sd: standard deviation of effect sizes 
:type sd: float
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

    # A simple case
    import fwdpy11
    # s N(0,0.1) and co-dominant
    gaussianS = fwdpy11.GaussianS(0,1,1,0.1,1)
)delim";

void
init_GaussianS(py::module& m)
{
    py::class_<fwdpy11::GaussianS, fwdpy11::Sregion>(
        m, "GaussianS", "Gaussian distribution of effect sizes")
        .def(py::init([](double beg, double end, double weight, double sd, double h,
                         bool coupled, std::uint16_t label, double scaling) {
                 return fwdpy11::GaussianS(
                     fwdpy11::Region(beg, end, weight, coupled, label), scaling, sd, h);
             }),
             py::arg("beg"), py::arg("end"), py::arg("weight"), py::arg("sd"),
             py::arg("h") = 1.0, py::arg("coupled") = true, py::arg("label") = 0,
             py::arg("scaling") = 1.0, INIT_DOCSTRING)
        .def_readonly("sd", &fwdpy11::GaussianS::sd)
        .def_readonly("h", &fwdpy11::GaussianS::dominance)
        .def("__repr__", &fwdpy11::GaussianS::repr)
        .def(py::pickle([](const fwdpy11::GaussianS& self) { return self.pickle(); },
                        [](py::tuple t) { return fwdpy11::GaussianS::unpickle(t); }));
}

