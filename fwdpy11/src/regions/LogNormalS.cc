#include <pybind11/pybind11.h>
#include <fwdpy11/regions/LogNormalS.hpp>

namespace py = pybind11;

static const auto CLASS_DOCSTRING = R"delim(
Log-normal distribution of effect sizes.

.. versionadded:: 0.7.0
)delim";

static const auto INIT_DOCSTRING =
    R"delim(
Constructor

:param beg: the beginning of the region
:type beg: float
:param end: the end of the region
:type end: float
:param weight: the weight to assign
:type weight: float
:param zeta: the zeta parameter
:type zeta: float
:param sigma: the sigma parameter
:type sigma: float
:param h: the dominance
:type h: float
:param coupled: if True, the weight is converted to(end-beg)*weight
:type Constructor: bool
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
    gdist=fwdpy11.LogNormalS(0, 1, 1, -0.1, 0.35)
)delim";

static const auto MV_DOCSTRING = R"delim(
Create an instance compatible with 
:class:`fwdpy11.mvDES`. See :ref:`mvdes` for
details.

:param beg: the beginning of the region
:type beg: float
:param end: the end of the region
:type end: float
:param weight: the weight to assign
:type weight: float
:param h: the dominance
:type h: float
:param coupled: if True, the weight is converted to(end-beg)*weight
:type coupled: bool
:param label: Fill :attr:`fwdpy11.Mutation.label` with this value.
:type label: numpy.uint16
:param scaling: The scaling of the DFE
:type scaling: float
)delim";

void
init_LogNormalS(py::module& m)
{
    py::class_<fwdpy11::LogNormalS, fwdpy11::Sregion>(m, "LogNormalS", CLASS_DOCSTRING)
        .def(
            py::init([](double beg, double end, double weight, double zeta, double sigma,
                        double h, bool coupled, std::uint16_t label, double scaling) {
                return fwdpy11::LogNormalS(
                    fwdpy11::Region(beg, end, weight, coupled, label), scaling, zeta,
                    sigma, h);
            }),
            py::arg("beg"), py::arg("end"), py::arg("weight"), py::arg("zeta"),
            py::arg("sigma"), py::arg("h") = 1.0, py::arg("coupled") = true,
            py::arg("label") = 0, py::arg("scaling") = 1.0, INIT_DOCSTRING)
        .def_static(
            "mv",
            [](double beg, double end, double weight, double h, bool coupled,
               std::uint16_t label, double scaling) {
                return fwdpy11 ::LogNormalS(
                    fwdpy11::Region(beg, end, weight, coupled, label), scaling, h);
            },
            py::arg("beg"), py::arg("end"), py::arg("weight"), py::arg("h") = 1.0,
            py::arg("coupled") = true, py::arg("label") = 0, py::arg("scaling") = 1.0,
            MV_DOCSTRING)
        .def_readonly("zeta", &fwdpy11::LogNormalS::zeta)
        .def_readonly("sigma", &fwdpy11::LogNormalS::sigma)
        .def_readonly("h", &fwdpy11::LogNormalS::dominance)
        .def("__repr__", &fwdpy11::LogNormalS::repr)
        .def(py::pickle([](const fwdpy11::LogNormalS& self) { return self.pickle(); },
                        [](py::tuple t) { return fwdpy11::LogNormalS::unpickle(t); }));
}
