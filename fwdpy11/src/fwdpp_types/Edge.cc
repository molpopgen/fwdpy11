#include <fwdpp/ts/edge.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <sstream>

namespace py = pybind11;

void
init_ts_Edge(py::module& m)
{
    PYBIND11_NUMPY_DTYPE(fwdpp::ts::edge, left, right, parent, child);
    py::class_<fwdpp::ts::edge>(m, "Edge",
                                R"delim(
            An edge in a tree sequence.  An edge
            represents the transmission of the
            half-open genomic interval [left,right) from
            parent to child.

            .. versionadded:: 0.2.0
            )delim")
        .def_readonly("left", &fwdpp::ts::edge::left,
                      "Left edge of interval, inclusive.")
        .def_readonly("right", &fwdpp::ts::edge::right,
                      "Right edge of interval, exclusive.")
        .def_readonly("parent", &fwdpp::ts::edge::parent, "Node id of parent")
        .def_readonly("child", &fwdpp::ts::edge::child, "Node id of child")
        .def("__repr__",
             [](const fwdpp::ts::edge& self) {
                 std::ostringstream out;
                 out.precision(4);
                 out << "Edge(parent=" << self.parent
                     << ", child=" << self.child << ", left=" << self.left
                     << ", right=" << self.right << ")";
                 return out.str();
             })
        .def(py::pickle(
            [](const fwdpp::ts::edge& e) {
                return py::make_tuple(e.left, e.right, e.parent, e.child);
            },
            [](py::tuple t) {
                return fwdpp::ts::edge{
                    t[0].cast<decltype(fwdpp::ts::edge::left)>(),
                    t[1].cast<decltype(fwdpp::ts::edge::right)>(),
                    t[2].cast<decltype(fwdpp::ts::edge::parent)>(),
                    t[3].cast<decltype(fwdpp::ts::edge::child)>()
                };
            }));
}

