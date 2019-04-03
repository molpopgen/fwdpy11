#include <pybind11/pybind11.h>
#include <fwdpy11/regions/UniformS.hpp>

namespace py = pybind11;

void
init_UniformS(py::module& m)
{
    py::class_<fwdpy11::UniformS, fwdpy11::Sregion>(
        m, "UniformS", "Uniform distrubution of effect sizes")
        .def(py::init([](double beg, double end, double weight, double lo,
                         double hi, double h, bool coupled,
                         std::uint16_t label, double scaling) {
                 return fwdpy11::UniformS(
                     fwdpy11::Region(beg, end, weight, coupled, label),
                     scaling, lo, hi, h);
             }),
             py::arg("beg"), py::arg("end"), py::arg("weight"), py::arg("lo"),
             py::arg("hi"), py::arg("h") = 1.0, py::arg("coupled") = true,
             py::arg("label") = 0, py::arg("scaling") = 1.0,
             R"delim(
            Constructor

            :param beg: the beginning of the region
            :param end: the end of the region
            :param weight: the weight to assign
            :param lo: lower bound on s
            :param hi: upper bound on s
            :param h: the dominance
            :param coupled: if True, the weight is converted to (end-beg)*weight
            :param label: Not relevant to recombining regions.
                Otherwise, this value will be used to take
                mutations from this region.
            :param scaling: The scaling of the DFE

            When coupled is True, the "weight" may be
            interpreted as a "per base pair"
            (or per unit, generally speaking) term.

            Example:

            .. testcode::

                #A simple case
                import fwdpy11
                #s is uniform on [-1, 0]
                uniformS = fwdpy11.UniformS(0,1,1,-1,0,0)
            )delim")
        .def_readonly("lo", &fwdpy11::UniformS::lo)
        .def_readonly("hi", &fwdpy11::UniformS::hi)
        .def_readonly("h", &fwdpy11::UniformS::dominance)
        .def("__repr__",&fwdpy11::UniformS::repr)
        .def(py::pickle(
            [](const fwdpy11::UniformS& self) { return self.pickle(); },
            [](py::tuple t) { return fwdpy11::UniformS::unpickle(t); }));
}

