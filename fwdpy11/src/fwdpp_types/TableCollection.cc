#include <fwdpp/ts/std_table_collection.hpp>
#include <fwdpy11/util/convert_lists.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(fwdpp::ts::std_table_collection::edge_table);
PYBIND11_MAKE_OPAQUE(fwdpp::ts::std_table_collection::node_table);
PYBIND11_MAKE_OPAQUE(fwdpp::ts::std_table_collection::mutation_table);
PYBIND11_MAKE_OPAQUE(fwdpp::ts::std_table_collection::site_table);

namespace
{
    template <typename T>
    void
    swap_with_empty(T& t)
    {
        T temp;
        t.swap(temp);
    }
}

void
init_ts_TableCollection(py::module& m)
{
    // A table_collection will not be user-constructible.  Rather,
    // they will be members of the Population classes.
    // TODO: work out conversion to msprime format
    // TODO: allow access to the "right" member functions
    py::class_<fwdpp::ts::std_table_collection,
               std::shared_ptr<fwdpp::ts::std_table_collection>>(m, "ll_TableCollection")
        .def(py::init([](std::shared_ptr<fwdpp::ts::std_table_collection>& self) {
            if (self == nullptr)
                {
                    throw std::invalid_argument("input tables cannot be nullptr");
                }
            return self;
        }))
        .def_readonly("_edges", &fwdpp::ts::std_table_collection::edges)
        .def_readonly("_nodes", &fwdpp::ts::std_table_collection::nodes)
        .def_readonly("_mutations", &fwdpp::ts::std_table_collection::mutations)
        .def_readonly("_sites", &fwdpp::ts::std_table_collection::sites)
        .def_readonly("_input_left", &fwdpp::ts::std_table_collection::input_left)
        .def_readonly("_output_right", &fwdpp::ts::std_table_collection::output_right)
        .def("_clear_edges",
             [](fwdpp::ts::std_table_collection& self) { swap_with_empty(self.edges); })
        .def("_clear_nodes",
             [](fwdpp::ts::std_table_collection& self) { swap_with_empty(self.nodes); })
        .def("_clear_sites",
             [](fwdpp::ts::std_table_collection& self) { swap_with_empty(self.sites); })
        .def("_clear_mutations",
             [](fwdpp::ts::std_table_collection& self) {
                 swap_with_empty(self.mutations);
             })
        .def("_clear_indexes",
             [](fwdpp::ts::std_table_collection& self) {
                 swap_with_empty(self.input_left);
                 swap_with_empty(self.output_right);
             })
        .def("_build_indexes", &fwdpp::ts::std_table_collection::build_indexes)
        .def_property_readonly("_genome_length",
                               &fwdpp::ts::std_table_collection::genome_length)
        .def("__eq__",
             [](const fwdpp::ts::std_table_collection& lhs,
                const fwdpp::ts::std_table_collection& rhs) { return lhs == rhs; })
        .def(py::pickle(
            [](const fwdpp::ts::std_table_collection& tables) {
                return py::make_tuple(tables.genome_length(),
                                      fwdpy11::vector_to_list(tables.nodes),
                                      fwdpy11::vector_to_list(tables.edges),
                                      fwdpy11::vector_to_list(tables.mutations),
                                      fwdpy11::vector_to_list(tables.sites));
            },
            [](py::tuple t) {
                auto length = t[0].cast<double>();
                fwdpp::ts::std_table_collection tables(length);
                tables.nodes = fwdpy11::list_to_vector<
                    fwdpp::ts::std_table_collection::node_table>(t[1].cast<py::list>());
                tables.edges = fwdpy11::list_to_vector<
                    fwdpp::ts::std_table_collection::edge_table>(t[2].cast<py::list>());
                tables.mutations = fwdpy11::list_to_vector<
                    fwdpp::ts::std_table_collection::mutation_table>(
                    t[3].cast<py::list>());
                tables.sites = fwdpy11::list_to_vector<
                    fwdpp::ts::std_table_collection::site_table>(t[4].cast<py::list>());
                tables.build_indexes();
                return tables;
            }));
}

