#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <fwdpp/ts/table_collection.hpp>
#include <fwdpp/ts/table_simplifier.hpp>

#include <fwdpy11/types/Population.hpp>
#include <fwdpy11/rng.hpp>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(fwdpp::ts::edge_vector);
PYBIND11_MAKE_OPAQUE(fwdpp::ts::node_vector);
PYBIND11_MAKE_OPAQUE(fwdpp::ts::mutation_key_vector);

PYBIND11_MODULE(ts, m)
{
    // TODO: how do I docstring this?
    m.attr("NULL_NODE") = py::int_(fwdpp::ts::TS_NULL_NODE);

    // The low-level types are numpy dtypes
    PYBIND11_NUMPY_DTYPE(fwdpp::ts::node, population, time);
    PYBIND11_NUMPY_DTYPE(fwdpp::ts::edge, left, right, parent, child);
    PYBIND11_NUMPY_DTYPE(fwdpp::ts::mutation_record, node, key);

    // The low level types are also Python classes
    py::class_<fwdpp::ts::node>(m, "Node")
        .def_readonly("population", &fwdpp::ts::node::population)
        .def_readonly("time", &fwdpp::ts::node::time);

    py::class_<fwdpp::ts::edge>(m, "Edge")
        .def_readonly("left", &fwdpp::ts::edge::left)
        .def_readonly("right", &fwdpp::ts::edge::right)
        .def_readonly("parent", &fwdpp::ts::edge::parent)
        .def_readonly("child", &fwdpp::ts::edge::child);

    py::class_<fwdpp::ts::mutation_record>(m, "MutationRecord")
        .def_readonly("node", &fwdpp::ts::mutation_record::node)
        .def_readonly("key", &fwdpp::ts::mutation_record::key);

    // indexed_edge cannot be a dtype.  That's probably ok,
    // as no-one will be using them for purposes other than viewing?
    // For now, I won't even bind the C++ vector of these...
    py::class_<fwdpp::ts::indexed_edge>(
        m, "IndexedEdge",
        "An edge keyed for efficient traversal of tree sequences.")
        .def_readonly("pos", &fwdpp::ts::indexed_edge::pos)
        .def_readonly("time", &fwdpp::ts::indexed_edge::time)
        .def_readonly("parent", &fwdpp::ts::indexed_edge::parent)
        .def_readonly("child", &fwdpp::ts::indexed_edge::child);

    // The tables are visible w/o copy via numpy
    py::bind_vector<fwdpp::ts::edge_vector>(
        m, "EdgeTable", py::buffer_protocol(), py::module_local(false));
    py::bind_vector<fwdpp::ts::node_vector>(
        m, "NodeTable", py::buffer_protocol(), py::module_local(false));
    py::bind_vector<fwdpp::ts::mutation_key_vector>(
        m, "MutationTable", py::buffer_protocol(), py::module_local(false));

    // A table_collection will not be user-constructible.  Rather,
    // they will be members of the Population classes.
    // TODO: work out conversion to msprime format
    // TODO: allow preserved_nodes to be cleared
    // TODO: allow access to the "right" member functions
    py::class_<fwdpp::ts::table_collection>(m, "TableCollection")
        .def_readonly("L", &fwdpp::ts::table_collection::L)
        .def_readonly("edges", &fwdpp::ts::table_collection::edge_table)
        .def_readonly("nodes", &fwdpp::ts::table_collection::node_table)
        .def_readonly("mutations",
                      &fwdpp::ts::table_collection::mutation_table)
        .def_readonly("input_left", &fwdpp::ts::table_collection::input_left)
        .def_readonly("output_right",
                      &fwdpp::ts::table_collection::output_right);

    py::class_<fwdpp::ts::table_simplifier>(m, "TableSimplifier")
        .def(py::init<double>(), py::arg("region_length"))
        .def(py::init([](const fwdpp::ts::table_collection& tables) {
            return fwdpp::ts::table_simplifier(tables.L);
        }))
        .def("__call__",
             [](fwdpp::ts::table_simplifier& simplifier,
                const fwdpy11::Population& pop,
                fwdpp::ts::table_collection& tables,
                const std::vector<fwdpp::ts::TS_NODE_INT>& samples) {
                 return simplifier.simplify(tables, samples, pop.mutations);
             })
        .def("simplify_copy",
             [](fwdpp::ts::table_simplifier& simplifier,
                const fwdpy11::Population& pop,
                const fwdpp::ts::table_collection& tables,
                const std::vector<fwdpp::ts::TS_NODE_INT>& samples) {
                 auto t(tables);
                 auto idmap = simplifier.simplify(t, samples, pop.mutations);
                 return py::make_tuple(std::move(t), std::move(idmap));
             });

    m.def("infinite_sites",
          [](const fwdpy11::GSLrng_t& rng, fwdpp::ts::table_collection& tables,
             const std::vector<fwdpp::ts::TS_NODE_INT>& samples, const double mu) {
              //TODO: implement
          });
}
