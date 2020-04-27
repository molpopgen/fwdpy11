#include <sstream>
#include <fwdpp/ts/indexed_edge.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

void
init_ts_IndexedEdge(py::module& m)
{
    // indexed_edge cannot be a dtype because it has a constructor.
    // That's probably ok, as no-one will be using them for purposes other than viewing?
    // For now, I won't even bind the C++ vector of these...
    py::class_<fwdpp::ts::indexed_edge>(
        m, "IndexedEdge",
        R"delim(An edge keyed for efficient traversal of tree sequences.
        
        .. versionadded:: 0.5.1
        )delim")
        .def_readonly("pos", &fwdpp::ts::indexed_edge::pos)
        .def_readonly("time", &fwdpp::ts::indexed_edge::time)
        .def_readonly("parent", &fwdpp::ts::indexed_edge::parent)
        .def_readonly("child", &fwdpp::ts::indexed_edge::child)
        .def("__repr__", [](const fwdpp::ts::indexed_edge& self) {
            std::ostringstream out;
            out.precision(4);
            out << "IndexedEdge(pos=" << self.pos << ", time=" << self.time
                << ", parent=" << self.parent << ", child=" << self.child
                << ')';
            return out.str();
        });
}

