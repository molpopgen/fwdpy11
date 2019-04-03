#include <pybind11/pybind11.h>
#include <fwdpy11/regions/ConstantS.hpp>

namespace py = pybind11;

void
init_ConstantS(py::module& m)
{
    py::class_<fwdpy11::ConstantS, fwdpy11::Sregion>(
        m, "ConstantS", "Mutations with fixed effect sizes")
        .def(py::init([](double beg, double end, double weight, double s,
                         double h, bool coupled, std::uint16_t label,
                         double scaling) {
                 return fwdpy11::ConstantS(
                     fwdpy11::Region(beg, end, weight, coupled, label),
                     scaling, s, h);
             }),
             py::arg("beg"), py::arg("end"), py::arg("weight"), py::arg("s"),
             py::arg("h") = 1.0, py::arg("coupled") = true,
             py::arg("label") = 0, py::arg("scaling") = 1.0,
             R"delim(
            Constructor

            :param beg: the beginning of the region
            :param end: the end of the region
            :param weight: the weight to assign
            :param s: the selection coefficient
            :param h: the dominance
            :param coupled: if True, the weight is converted to (end-beg)*weight
            :param label: Not relevant to recombining regions.
                Otherwise, this value will be used to take mutations
                from this region.
            :param scaling: The scaling of the DFE

            When coupled is True, the "weight" may be interpreted
            as a "per base pair"
            (or per unit, generally speaking) term.

            Example:

            .. testcode::

                #A simple case
                import fwdpy11
                #s = -0.1 and h = 0
                constantS = fwdpy11.ConstantS(0,1,1,-0.1,0)
            )delim")
        .def_readonly("esize", &fwdpy11::ConstantS::esize)
        .def_readonly("s", &fwdpy11::ConstantS::esize)
        .def_readonly("h", &fwdpy11::ConstantS::dominance)
        .def("__repr__",&fwdpy11::ConstantS::repr)
        .def(py::pickle(
            [](const fwdpy11::ConstantS& self) { return self.pickle(); },
            [](py::tuple t) { return fwdpy11::ConstantS::unpickle(t); }));
}

