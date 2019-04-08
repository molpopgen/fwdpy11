#include <fwdpp/ts/node.hpp>
#include <fwdpy11/util/convert_lists.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<fwdpp::ts::node>);

void
init_ts_NodeTable(py::module& m)
{
    py::bind_vector<std::vector<fwdpp::ts::node>>(
        m, "NodeTable", py::buffer_protocol(), py::module_local(false),
        R"delim(
        An NodeTable is a container of :class:`fwdpy11.Node`.

        An NodeTable supports the Python buffer protocol, allowing data
        to be viewed in Python as a numpy record array.  No copy is made
        when generating such views.

        .. versionadded:: 0.2.0
        )delim")
        .def(py::pickle(
            [](const std::vector<fwdpp::ts::node>& nodes) {
                return fwdpy11::vector_to_list(nodes);
            },
            [](py::list l) {
                return fwdpy11::list_to_vector<std::vector<fwdpp::ts::node>>(
                    l);
            }));
}

