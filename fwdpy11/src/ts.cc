#include <algorithm>
#include <vector>
#include <stdexcept>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <fwdpp/ts/table_collection.hpp>
#include <fwdpp/ts/marginal_tree.hpp>
#include <fwdpp/ts/tree_visitor.hpp>
#include <fwdpp/ts/recycling.hpp>
#include <fwdpp/ts/mutate_tables.hpp>
#include <fwdpp/ts/count_mutations.hpp>

#include <fwdpy11/types/Population.hpp>
#include <fwdpy11/rng.hpp>
#include <fwdpy11/policies/mutation.hpp>
#include <fwdpy11/numpy/array.hpp>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(fwdpp::ts::edge_vector);
PYBIND11_MAKE_OPAQUE(fwdpp::ts::node_vector);
PYBIND11_MAKE_OPAQUE(fwdpp::ts::mutation_key_vector);

void init_ts(py::module&);

inline std::size_t
generate_neutral_variants(fwdpp::flagged_mutation_queue& recycling_bin,
                          fwdpy11::Population& pop,
                          const fwdpy11::GSLrng_t& rng, const double left,
                          const double right, const fwdpp::uint_t generation)
{
    const auto uniform = [left, right, &rng]() {
        return gsl_ran_flat(rng.get(), left, right);
    };
    const auto return_zero = []() { return 0.0; };
    return fwdpy11::infsites_Mutation(recycling_bin, pop.mutations,
                                      pop.mut_lookup, generation, uniform,
                                      return_zero, return_zero, 0);
}

template <typename T>
py::list
vector_to_list(const T& t)
{
    py::list rv;
    for (auto& i : t)
        {
            rv.append(i);
        }
    return rv;
}

template <typename T>
T
list_to_vector(py::list l)
{
    T rv;
    rv.reserve(l.size());
    for (auto& i : l)
        {
            rv.push_back(i.cast<typename T::value_type>());
        }
    return rv;
}

PYBIND11_MODULE(ts, m)
{

    // The low-level types are numpy dtypes


    // The tables are visible w/o copy via numpy




    py::class_<fwdpp::ts::marginal_tree>(
        m, "MarginalTree",
        R"delim(A sparse tree representation of a non-recombining genomic segment.
        See :ref:`ts_data_types` for details.
        
        .. deprecated:: 0.3.0

            Deprecated in favor of :class:`fwdpy11.ts.TreeIterator`
        )delim")
        .def_readonly("left", &fwdpp::ts::marginal_tree::left,
                      "Left edge of genomic interval (inclusive)")
        .def_readonly("right", &fwdpp::ts::marginal_tree::right,
                      "Right edge of genomic interval (exclusive")
        .def_property_readonly("parents",
                               [](const fwdpp::ts::marginal_tree& m) {
                                   return fwdpy11::make_1d_ndarray(m.parents);
                               },
                               "Vector of child -> parent relationships")
        .def_property_readonly("leaf_counts",
                               [](const fwdpp::ts::marginal_tree& m) {
                                   return fwdpy11::make_1d_ndarray(
                                       m.leaf_counts);
                               },
                               "Leaf counts for each node")
        .def_property_readonly("preserved_leaf_counts",
                               [](const fwdpp::ts::marginal_tree& m) {
                                   return fwdpy11::make_1d_ndarray(
                                       m.preserved_leaf_counts);
                               },
                               "Ancient sample leaf counts for each node")
        .def_property_readonly("left_sib",
                               [](const fwdpp::ts::marginal_tree& m) {
                                   return fwdpy11::make_1d_ndarray(m.left_sib);
                               },
                               "Return the left sibling of the current node")
        .def_property_readonly("right_sib",
                               [](const fwdpp::ts::marginal_tree& m) {
                                   return fwdpy11::make_1d_ndarray(
                                       m.right_sib);
                               },
                               "Return the right sibling of the current node")
        .def_property_readonly("left_child",
                               [](const fwdpp::ts::marginal_tree& m) {
                                   return fwdpy11::make_1d_ndarray(
                                       m.left_child);
                               },
                               "Mapping of current node id to its left child")
        .def_property_readonly("right_child",
                               [](const fwdpp::ts::marginal_tree& m) {
                                   return fwdpy11::make_1d_ndarray(
                                       m.right_child);
                               },
                               "Mapping of current node id to its right child")
        .def_property_readonly("left_sample",
                               [](const fwdpp::ts::marginal_tree& m) {
                                   return fwdpy11::make_1d_ndarray(
                                       m.left_sample);
                               })
        .def_property_readonly("right_sample",
                               [](const fwdpp::ts::marginal_tree& m) {
                                   return fwdpy11::make_1d_ndarray(
                                       m.right_sample);
                               })
        .def_property_readonly("next_sample",
                               [](const fwdpp::ts::marginal_tree& m) {
                                   return fwdpy11::make_1d_ndarray(
                                       m.next_sample);
                               })
        .def_property_readonly("sample_index_map",
                               [](const fwdpp::ts::marginal_tree& m) {
                                   return fwdpy11::make_1d_ndarray(
                                       m.sample_index_map);
                               })
        .def_readonly("sample_size", &fwdpp::ts::marginal_tree::sample_size)
        .def("total_time",
             [](const fwdpp::ts::marginal_tree& m,
                const fwdpp::ts::node_vector& nodes) {
                 if (m.parents.size() != nodes.size())
                     {
                         throw std::invalid_argument(
                             "node table length does not equal number of "
                             "nodes in marginal tree");
                     }
                 double tt = 0.0;
                 for (std::size_t i = 0; i < m.parents.size(); ++i)
                     {
                         if (m.parents[i] != fwdpp::ts::TS_NULL_NODE)
                             {
                                 tt += nodes[i].time
                                       - nodes[m.parents[i]].time;
                             }
                     }
                 return tt;
             },
             "Return the sum of branch lengths");

    py::class_<fwdpp::ts::tree_visitor>(
        m, "TreeVisitor",
        R"delim(Allows left-to-right visiting of the marginal trees embedded in a
        :class:`fwdpy11.ts.TableCollection`
        
        .. deprecated:: 0.3.0

            Deprecated in favor of :class:`fwdpy11.ts.TreeIterator`
        )delim")
        .def(py::init<const fwdpp::ts::table_collection&,
                      const std::vector<fwdpp::ts::TS_NODE_INT>&>(),
             py::arg("tables"), py::arg("samples"))
        .def(py::init<const fwdpp::ts::table_collection&,
                      const std::vector<fwdpp::ts::TS_NODE_INT>&,
                      const std::vector<fwdpp::ts::TS_NODE_INT>&>(),
             py::arg("tables"), py::arg("samples"), py::arg("ancient_samples"))
        .def("tree", &fwdpp::ts::tree_visitor::tree,
             py::return_value_policy::reference_internal,
             "Obtain the current marginal tree as an instance of "
             ":class:`fwdpy11.ts.MarginalTree")
        .def("__call__",
             [](fwdpp::ts::tree_visitor& tv, const bool update_samples) {
                 bool rv = false;
                 if (update_samples)
                     {
                         rv = tv(std::true_type(), std::true_type());
                     }
                 else
                     {
                         rv = tv(std::true_type(), std::false_type());
                     }
                 return rv;
             },
             py::arg("update_samples"),
             R"delim(
             Advance to the next tree.

             :param update_samples: If True, update the "samples list"
             :type update_samples: booliean

             By default, leaf counts and ancestral leaf counts are
             always updated. The latter is only updated if you separated
             modern from ancient samples when constructing the object.
             )delim");

    m.def("infinite_sites", [](const fwdpy11::GSLrng_t& rng,
                               fwdpy11::Population& pop, const double mu) {
        fwdpp::flagged_mutation_queue recycling_bin
            = fwdpp::ts::make_mut_queue(pop.mcounts,
                                        pop.mcounts_from_preserved_nodes);
        std::vector<fwdpp::ts::TS_NODE_INT> samples(2 * pop.N);
        std::iota(samples.begin(), samples.end(), 0);
        const auto apply_mutations =
            [&recycling_bin, &rng, &pop](const double left, const double right,
                                         const fwdpp::uint_t generation) {
                return generate_neutral_variants(recycling_bin, pop, rng, left,
                                                 right, generation);
            };
        auto nmuts = fwdpp::ts::mutate_tables(rng, apply_mutations, pop.tables,
                                              samples, mu);
        std::sort(
            pop.tables.mutation_table.begin(), pop.tables.mutation_table.end(),
            [&pop](const fwdpp::ts::mutation_record& a,
                   const fwdpp::ts::mutation_record& b) {
                return pop.mutations[a.key].pos < pop.mutations[b.key].pos;
            });
        fwdpp::ts::count_mutations(pop.tables, pop.mutations, samples,
                                   pop.mcounts,
                                   pop.mcounts_from_preserved_nodes);
        return nmuts;
    });

    init_ts(m);
}
