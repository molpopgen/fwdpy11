#pragma GCC visibility push(default)
#include <fwdpp/ts/exceptions.hpp>
#pragma GCC visibility pop
#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_mutation_base(py::module &);
void init_HaploidGenome(py::module &);
void init_data_matrix(py::module &);
void init_ts_Node(py::module &);
void init_ts_Edge(py::module &);
void init_ts_MutationRecord(py::module &);
void init_ts_Site(py::module &m);
void init_NULL_NODE(py::module &);
void init_ts_NodeTable(py::module &m);
void init_ts_EdgeTable(py::module &m);
void init_ts_MutationTable(py::module &);
void init_ts_SiteTable(py::module &);
void init_ts_TableCollection(py::module &);

void
initialize_fwdpp_types(py::module &m)
{
    py::register_exception<fwdpp::ts::tables_error>(m, "TablesError");
    py::register_exception<fwdpp::ts::samples_error>(m, "SamplesError");


    // For testing only, mostly to check symbol visibility of 
    // exception types, which has been an issue on macOS
    m.def("_throw_TablesError",
          []() { throw fwdpp::ts::tables_error("this is a TablesError"); });
    m.def("_throw_SamplesError",
          []() { throw fwdpp::ts::samples_error("this is a SamplesError"); });

    init_mutation_base(m);
    init_HaploidGenome(m);
    init_data_matrix(m);

    // Types related to tree sequenc recording
    init_NULL_NODE(m);
    init_ts_Node(m);
    init_ts_Edge(m);
    init_ts_MutationRecord(m);
    init_ts_Site(m);
    init_ts_NodeTable(m);
    init_ts_EdgeTable(m);
    init_ts_MutationTable(m);
    init_ts_SiteTable(m);
    init_ts_TableCollection(m);
}
