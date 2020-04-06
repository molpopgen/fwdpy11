#include <pybind11/pybind11.h>
#include <fwdpy11/regions/LogNormalS.hpp>

namespace py = pybind11;

void
init_LogNormalS(py::module& m)
{
    py::class_<fwdpy11::LogNormalS, fwdpy11::Sregion>(
        m, "LogNormalS", "Log-normal distribution of effect sizes")
        .def(
            py::init([](double beg, double end, double weight, double zeta, double sigma,
                        double h, bool coupled, std::uint16_t label, double scaling) {
                return fwdpy11::LogNormalS(
                    fwdpy11::Region(beg, end, weight, coupled, label), scaling, zeta,
                    sigma, h);
            }),
            py::arg("beg"), py::arg("end"), py::arg("weight"), py::arg("zeta"),
            py::arg("sigma"), py::arg("h") = 1.0, py::arg("coupled") = true,
            py::arg("label") = 0, py::arg("scaling") = 1.0,
            R"delim(
            Constructor

            :param beg: the beginning of the region
            :param end: the end of the region
            :param weight: the weight to assign
            :param zeta: the zeta parameter
            :param sigma: the sigma parameter
            :param h: the dominance
            :param coupled: if True, the weight is converted to(end-beg)*weight
            :param label: Not relevant to recombining regions.
                Otherwise, this value will be used
                to take mutations from this region.
            :param scaling: The scaling of the DFE

            When coupled is True, the "weight" may be
            interpreted as a "per base pair"
            (or per unit, generally speaking) term.

            Example:

            .. testcode::

                # A simple case
                import fwdpy11
                gdist=fwdpy11.LogNormalS(0, 1, 1, -0.1, 0.35)

            .. versionadded: : 0.7.0
        )delim")
        .def_static(
            "mv",
            [](double beg, double end, double weight, double h, bool coupled,
               std::uint16_t label, double scaling) {
                return fwdpy11 ::LogNormalS(
                    fwdpy11::Region(beg, end, weight, coupled, label), scaling, h);
            },
            py::arg("beg"), py::arg("end"), py::arg("weight"), py::arg("h") = 1.0,
            py::arg("coupled") = true, py::arg("label") = 0, py::arg("scaling") = 1.0)
        .def_readonly("zeta", &fwdpy11::LogNormalS::zeta)
        .def_readonly("sigma", &fwdpy11::LogNormalS::sigma)
        .def_readonly("h", &fwdpy11::LogNormalS::dominance)
        .def("__repr__", &fwdpy11::LogNormalS::repr)
        .def(py::pickle([](const fwdpy11::LogNormalS& self) { return self.pickle(); },
                        [](py::tuple t) { return fwdpy11::LogNormalS::unpickle(t); }));
}
