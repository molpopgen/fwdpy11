#include <fwdpp/ts/node.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <sstream>

namespace py = pybind11;

void
init_ts_Node(py::module& m)
{
    PYBIND11_NUMPY_DTYPE(fwdpp::ts::node, deme, time);

    // The low level types are also Python classes
    py::class_<fwdpp::ts::node>(m, "Node",
                                R"delim(
            A node in a tree sequence.

            .. versionadded:: 0.2.0
            )delim")
        .def_readonly("deme", &fwdpp::ts::node::deme,
                      R"delim(
            For models of discrete population structure,
            this field is the deme of the node.
            )delim")
        .def_readonly("time", &fwdpp::ts::node::time,
                      "Birth time of the node, recorded forwards in time.")
        .def("__repr__",
             [](const fwdpp::ts::node& self) {
                 std::ostringstream out;
                 out.precision(4);
                 out << "Node(time=" << self.time
                     << ", deme=" << self.deme << ")";
                 return out.str();
             })
        .def(py::pickle(
            [](const fwdpp::ts::node& n) {
                return py::make_tuple(n.deme, n.time);
            },
            [](py::tuple t) {
                return fwdpp::ts::node{
                    t[0].cast<decltype(fwdpp::ts::node::deme)>(),
                    t[1].cast<double>()
                };
            }));
}

