#include <fwdpp/ts/table_collection.hpp>
#include <fwdpy11/util/convert_lists.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(fwdpp::ts::edge_vector);
PYBIND11_MAKE_OPAQUE(fwdpp::ts::node_vector);
PYBIND11_MAKE_OPAQUE(fwdpp::ts::mutation_key_vector);
PYBIND11_MAKE_OPAQUE(fwdpp::ts::site_vector);

void
init_ts_TableCollection(py::module& m)
{
    // A table_collection will not be user-constructible.  Rather,
    // they will be members of the Population classes.
    // TODO: work out conversion to msprime format
    // TODO: allow preserved_nodes to be cleared
    // TODO: allow access to the "right" member functions
    py::class_<fwdpp::ts::table_collection>(
        m, "TableCollection",
        "A table collection representing a succinct tree sequence.")
        .def_property_readonly(
            "L", &fwdpp::ts::table_collection::genome_length, "Genome length")
        .def_readonly("edges", &fwdpp::ts::table_collection::edge_table,
                      "The :class:`fwdpy11.EdgeTable`.")
        .def_readonly("nodes", &fwdpp::ts::table_collection::node_table,
                      "The :class:`fwdpy11.NodeTable`.")
        .def_readonly("mutations",
                      &fwdpp::ts::table_collection::mutation_table,
                      "The :class:`fwdpy11.MutationTable`.")
        .def_readonly("sites",
                      &fwdpp::ts::table_collection::site_table,
                      "The :class:`fwdpy11.SiteTable`.")
        .def_readonly("input_left", &fwdpp::ts::table_collection::input_left)
        .def_readonly("output_right",
                      &fwdpp::ts::table_collection::output_right)
        .def_readonly("preserved_nodes",
                      &fwdpp::ts::table_collection::preserved_nodes,
                      "List of nodes corresponding to ancient samples.")
        .def_property_readonly("genome_length",
                               &fwdpp::ts::table_collection::genome_length,
                               "Return the genome/sequence length.")
        .def("__eq__",
             [](const fwdpp::ts::table_collection& lhs,
                const fwdpp::ts::table_collection& rhs) { return lhs == rhs; })
        .def(py::pickle(
            [](const fwdpp::ts::table_collection& tables) {
                return py::make_tuple(
                    tables.genome_length(),
                    fwdpy11::vector_to_list(tables.node_table),
                    fwdpy11::vector_to_list(tables.edge_table),
                    fwdpy11::vector_to_list(tables.mutation_table),
                    fwdpy11::vector_to_list(tables.site_table),
                    fwdpy11::vector_to_list(tables.preserved_nodes));
            },
            [](py::tuple t) {
                auto length = t[0].cast<double>();
                fwdpp::ts::table_collection tables(length);
                tables.node_table
                    = fwdpy11::list_to_vector<fwdpp::ts::node_vector>(
                        t[1].cast<py::list>());
                tables.edge_table
                    = fwdpy11::list_to_vector<fwdpp::ts::edge_vector>(
                        t[2].cast<py::list>());
                tables.mutation_table
                    = fwdpy11::list_to_vector<fwdpp::ts::mutation_key_vector>(
                        t[3].cast<py::list>());
                tables.site_table
                    = fwdpy11::list_to_vector<fwdpp::ts::site_vector>(
                        t[4].cast<py::list>());
                tables.preserved_nodes = fwdpy11::list_to_vector<decltype(
                    fwdpp::ts::table_collection::preserved_nodes)>(
                    t[5].cast<py::list>());
                tables.build_indexes();
                return tables;
            }));
}

