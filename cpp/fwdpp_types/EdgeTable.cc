#include <fwdpp/ts/edge.hpp>
#include <fwdpy11/util/convert_lists.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<fwdpp::ts::edge>);

void
init_ts_EdgeTable(py::module& m)
{
    py::bind_vector<std::vector<fwdpp::ts::edge>>(
        m, "EdgeTable", py::buffer_protocol(), py::module_local(false),
        R"delim(
        An EdgeTable is a container of :class:`fwdpy11.Edge`.

        An EdgeTable supports the Python buffer protocol, allowing data
        to be viewed in Python as a numpy record array.  No copy is made
        when generating such views.

        .. versionadded:: 0.2.0
        )delim")
        .def(py::pickle(
            [](const std::vector<fwdpp::ts::edge>& edges) {
                return fwdpy11::vector_to_list(edges);
            },
            [](py::list l) {
                return fwdpy11::list_to_vector<std::vector<fwdpp::ts::edge>>(l);
            }));
}

