#include <fwdpp/ts/mutation_record.hpp>
#include <fwdpy11/util/convert_lists.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<fwdpp::ts::mutation_record>);

void
init_ts_MutationTable(py::module& m)
{
    py::bind_vector<std::vector<fwdpp::ts::mutation_record>>(
        m, "MutationTable", py::buffer_protocol(), py::module_local(false),
        R"delim(
        An MutationTable is a container of :class:`fwdpy11.MutationRecord`.

        An MutationTable supports the Python buffer protocol, allowing data
        to be viewed in Python as a numpy record array.  No copy is made
        when generating such views.

        .. versionadded:: 0.2.0
        )delim")
        .def(py::pickle(
            [](const std::vector<fwdpp::ts::mutation_record>& mutations) {
                return fwdpy11::vector_to_list(mutations);
            },
            [](py::list l) {
                return fwdpy11::list_to_vector<
                    std::vector<fwdpp::ts::mutation_record>>(l);
            }));
}

