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
void init_GeneticMapUnit(py::module &);
void init_PoissonInterval(py::module &);
void init_BinomialPoint(py::module &);
void init_PoissonPoint(py::module &);
void init_FixedCrossovers(py::module &);

void
initialize_fwdpp_types(py::module &m)
{
    init_mutation_base(m);
    init_HaploidGenome(m);
    init_data_matrix(m);

    init_GeneticMapUnit(m);
    init_PoissonInterval(m);
    init_BinomialPoint(m);
    init_PoissonPoint(m);
    init_FixedCrossovers(m);

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
