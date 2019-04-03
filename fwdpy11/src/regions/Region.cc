#include <pybind11/pybind11.h>
#include <fwdpy11/regions/Region.hpp>

namespace py = pybind11;

void
init_Region(py::module& m)
{
    py::class_<fwdpy11::Region>(
        m, "Region",
        "A genomic region, defined by half-open interval [beg, end)")
        .def(py::init<double, double, double, bool, std::uint16_t>(),
             py::arg("beg"), py::arg("end"), py::arg("weight"),
             py::arg("coupled") = true, py::arg("label") = 0,
             R"delim(
        Constructor

        :param beg: the beginning of the region
        :param end: the end of the region
        :param weight: the weight to assign
        :param coupled: if True, the weight is converted to (end-beg)*weight
        :param label: Not relevant to recombining regions.
            Otherwise, this value will be used to take mutations
            from this region.

        When coupled is True, the "weight" may be interpreted
        as a "per base pair"
        (or per unit, generally speaking) term.

        .. versionchanged:: 0.3.0

            Refactored from a pure Python class to a C++/pybind11 class

        Example:

        .. testcode::

            #A simple case
            import fwdpy11
            r = fwdpy11.Region(0,1,1)
            #A more "biological" case:
            #  The region covers positions 1 through 1,000,
            #  and the per-base pair "weight" is 1e-5:
            r = fwdpy11.Region(1,1000,1e-5,True)
        )delim")
        .def_readonly("b", &fwdpy11::Region::beg, "Beginning of region")
        .def_readonly("e", &fwdpy11::Region::end, "End of region")
        .def_readonly("w", &fwdpy11::Region::weight,
                      "Weight associated with the region")
        .def_readonly("c", &fwdpy11::Region::coupled,
                      "If True, total weight is a function of end-beg")
        .def_readonly("l", &fwdpy11::Region::label,
                      "Label associated with the region")
        .def("__repr__", &fwdpy11::Region::repr)
        .def(py::pickle(
            [](const fwdpy11::Region& self) { return self.pickle(); },
            [](py::tuple t) { return fwdpy11::Region::unpickle(t); }));
}

