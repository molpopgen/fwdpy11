#include <algorithm>
#include <vector>
#include <stdexcept>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <fwdpp/ts/table_collection.hpp>
#include <fwdpp/ts/table_simplifier.hpp>
#include <fwdpp/ts/marginal_tree.hpp>
#include <fwdpp/ts/tree_visitor.hpp>
#include <fwdpp/ts/recycling.hpp>
#include <fwdpp/ts/mutate_tables.hpp>
#include <fwdpp/ts/count_mutations.hpp>
#include <fwdpp/ts/generate_data_matrix.hpp>

#include <fwdpy11/types/Population.hpp>
#include <fwdpy11/rng.hpp>
#include <fwdpy11/policies/mutation.hpp>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(fwdpp::ts::edge_vector);
PYBIND11_MAKE_OPAQUE(fwdpp::ts::node_vector);
PYBIND11_MAKE_OPAQUE(fwdpp::ts::mutation_key_vector);

// TODO: need to clear out ancient nodes after
// copy.  The user should include them in "samples"
// if they are wanted.
py::tuple
simplify(const fwdpy11::Population& pop,
         const std::vector<fwdpp::ts::TS_NODE_INT>& samples)
{
    if (pop.tables.genome_length() == std::numeric_limits<double>::max())
        {
            throw std::invalid_argument(
                "population is not using tree sequences");
        }
    if (pop.tables.num_nodes() == 0)
        {
            throw std::invalid_argument(
                "population has empty TableCollection");
        }
    if (samples.empty())
        {
            throw std::invalid_argument("empty sample list");
        }
    if (std::any_of(samples.begin(), samples.end(),
                    [&pop](const fwdpp::ts::TS_NODE_INT s) {
                        return s == fwdpp::ts::TS_NULL_NODE
                               || static_cast<std::size_t>(s)
                                      >= pop.tables.num_nodes();
                    }))
        {
            throw std::invalid_argument("invalid sample list");
        }
    auto t(pop.tables);
    fwdpp::ts::table_simplifier simplifier(pop.tables.genome_length());
    auto rv = simplifier.simplify(t, samples, pop.mutations);
    t.build_indexes();
    return py::make_tuple(std::move(t), std::move(rv.first));
}

inline std::size_t
generate_neutral_variants(std::queue<std::size_t>& recycling_bin,
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
inline py::array_t<T>
make_1d_ndarray(const std::vector<T>& v)
{
    auto rv
        = py::array_t<T>({ v.size() }, { sizeof(T) }, v.data(), py::cast(v));
    return rv;
}

struct VariantIterator
{
    std::vector<fwdpp::ts::mutation_record>::const_iterator mbeg, mend;
    std::vector<double> pos;
    fwdpp::ts::tree_visitor tv;
    std::vector<std::int8_t> genotypes;
    VariantIterator(const fwdpp::ts::table_collection& tc,
                    const std::vector<fwdpy11::Mutation>& mutations,
                    const std::vector<fwdpp::ts::TS_NODE_INT>& samples)
        : mbeg(tc.mutation_table.begin()), mend(tc.mutation_table.end()),
          pos(), tv(tc, samples), genotypes(samples.size(), 0)
    {
        // Advance to first tree
        auto flag = tv(std::true_type(), std::true_type());
        if (flag == false)
            {
                throw std::invalid_argument(
                    "TableCollection contains no trees");
            }
        for (auto& m : mutations)
            {
                pos.push_back(m.pos);
            }
    }

    py::object
    next_variant()
    {
        if (!(mbeg < mend))
            {
                return py::none();
                throw py::stop_iteration();
            }
        auto m = tv.tree();
        if (pos[mbeg->key] >= m.right)
            {
                auto flag = tv(std::true_type(), std::false_type());
                if (flag == false)
                    {
                        return py::none();
                    }
            }
        std::fill(genotypes.begin(), genotypes.end(), 0);
        auto ls = m.left_sample[mbeg->node];
        if (ls != fwdpp::ts::TS_NULL_NODE)
            {
                auto rs = m.right_sample[mbeg->node];
                int nsteps = 1;
                while (true)
                    {
                        if (genotypes[ls] == 1)
                            {
                                throw std::runtime_error(
                                    "VariantIterator error");
                            }
                        genotypes[ls] = 1;
                        if (ls == rs)
                            {
                                break;
                            }
                        ls = m.next_sample[ls];
                        ++nsteps;
                    }
                if (nsteps != m.leaf_counts[mbeg->node])
                    {
                        throw std::runtime_error(
                            "VariantIterator: sample traversal error");
                    }
            }
        //py::print(std::distance(mbeg, mend), pos[mbeg->key], m.left, m.right);
        ++mbeg;
        return make_1d_ndarray(genotypes);
    }
};

PYBIND11_MODULE(ts, m)
{
    // TODO: how do I docstring this?
    m.attr("NULL_NODE") = py::int_(fwdpp::ts::TS_NULL_NODE);

    // The low-level types are numpy dtypes
    PYBIND11_NUMPY_DTYPE(fwdpp::ts::node, population, time);
    PYBIND11_NUMPY_DTYPE(fwdpp::ts::edge, left, right, parent, child);
    PYBIND11_NUMPY_DTYPE(fwdpp::ts::mutation_record, node, key);

    // The low level types are also Python classes
    py::class_<fwdpp::ts::node>(m, "Node",
                                R"delim(
            A node in a tree sequence.

            .. versionadded:: 0.2.0
            )delim")
        .def_readonly("population", &fwdpp::ts::node::population,
                      R"delim(
            For models of discrete population structure,
            this field is the population of the node.
            )delim")
        .def_readonly("time", &fwdpp::ts::node::time,
                      "Birth time of the node, recorded forwards in time.");

    py::class_<fwdpp::ts::edge>(m, "Edge",
                                R"delim(
            An edge in a tree sequence.  An edge
            represents the transmission of the
            half-open genomic interval [left,right) from
            parent to child.

            .. versionadded:: 0.2.0
            )delim")
        .def_readonly("left", &fwdpp::ts::edge::left,
                      "Left edge of interval, inclusive.")
        .def_readonly("right", &fwdpp::ts::edge::right,
                      "Right edge of interval, exclusive.")
        .def_readonly("parent", &fwdpp::ts::edge::parent, "Node id of parent")
        .def_readonly("child", &fwdpp::ts::edge::child, "Node id of child");

    py::class_<fwdpp::ts::mutation_record>(m, "MutationRecord",
                                           R"delim(
            A mutation entry in a tree sequence. A MutationRecord
            stores the node ID of the mutation and the index
            ("key") of the mutation in the population.

            .. versionadded:: 0.2.0
            )delim")
        .def_readonly("node", &fwdpp::ts::mutation_record::node,
                      "Node id of the mutation")
        .def_readonly("key", &fwdpp::ts::mutation_record::key,
                      "Index of the mutation in the population");

    // indexed_edge cannot be a dtype because it has a constructor.
    // That's probably ok, as no-one will be using them for purposes other than viewing?
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
        m, "EdgeTable", py::buffer_protocol(), py::module_local(false),
        R"delim(
        An EdgeTable is a container of :class:`fwdpy11.ts.Edge`.

        An EdgeTable supports the Python buffer protocol, allowing data
        to be viewed in Python as a numpy record array.  No copy is made
        when generating such views.

        .. versionadded:: 0.2.0
        )delim");

    py::bind_vector<fwdpp::ts::node_vector>(
        m, "NodeTable", py::buffer_protocol(), py::module_local(false),
        R"delim(
        An NodeTable is a container of :class:`fwdpy11.ts.Node`.

        An NodeTable supports the Python buffer protocol, allowing data
        to be viewed in Python as a numpy record array.  No copy is made
        when generating such views.

        .. versionadded:: 0.2.0
        )delim");

    py::bind_vector<fwdpp::ts::mutation_key_vector>(
        m, "MutationTable", py::buffer_protocol(), py::module_local(false),
        R"delim(
        An MutationTable is a container of :class:`fwdpy11.ts.MutationRecord`.

        An MutationTable supports the Python buffer protocol, allowing data
        to be viewed in Python as a numpy record array.  No copy is made
        when generating such views.

        .. versionadded:: 0.2.0
        )delim");

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
                      "The :class:`fwdpy11.ts.EdgeTable`.")
        .def_readonly("nodes", &fwdpp::ts::table_collection::node_table,
                      "The :class:`fwdpy11.ts.NodeTable`.")
        .def_readonly("mutations",
                      &fwdpp::ts::table_collection::mutation_table,
                      "The :class:`fwdpy11.ts.MutationTable`.")
        .def_readonly("input_left", &fwdpp::ts::table_collection::input_left)
        .def_readonly("output_right",
                      &fwdpp::ts::table_collection::output_right)
        .def_readonly("preserved_nodes",
                      &fwdpp::ts::table_collection::preserved_nodes,
                      "List of nodes corresponding to ancient samples.")
        .def("genome_length", &fwdpp::ts::table_collection::genome_length,
             "Return the genome/sequence length.");

    py::class_<fwdpp::ts::marginal_tree>(
        m, "MarginalTree",
        "A sparse tree representation of a non-recombining genomic segment. "
        "See :ref:`ts_data_types` for details.")
        .def_readonly("left", &fwdpp::ts::marginal_tree::left,
                      "Left edge of genomic interval (inclusive)")
        .def_readonly("right", &fwdpp::ts::marginal_tree::right,
                      "Right edge of genomic interval (exclusive")
        .def_property_readonly("parents",
                               [](const fwdpp::ts::marginal_tree& m) {
                                   return make_1d_ndarray(m.parents);
                               },
                               "Vector of child -> parent relationships")
        .def_property_readonly("leaf_counts",
                               [](const fwdpp::ts::marginal_tree& m) {
                                   return make_1d_ndarray(m.leaf_counts);
                               },
                               "Leaf counts for each node")
        .def_property_readonly("preserved_leaf_counts",
                               [](const fwdpp::ts::marginal_tree& m) {
                                   return make_1d_ndarray(
                                       m.preserved_leaf_counts);
                               },
                               "Ancient sample leaf counts for each node")
        .def_property_readonly("left_sib",
                               [](const fwdpp::ts::marginal_tree& m) {
                                   return make_1d_ndarray(m.left_sib);
                               },
                               "Return the left sibling of the current node")
        .def_property_readonly("right_sib",
                               [](const fwdpp::ts::marginal_tree& m) {
                                   return make_1d_ndarray(m.right_sib);
                               },
                               "Return the right sibling of the current node")
        .def_property_readonly("left_child",
                               [](const fwdpp::ts::marginal_tree& m) {
                                   return make_1d_ndarray(m.left_child);
                               },
                               "Mapping of current node id to its left child")
        .def_property_readonly("right_child",
                               [](const fwdpp::ts::marginal_tree& m) {
                                   return make_1d_ndarray(m.right_child);
                               },
                               "Mapping of current node id to its right child")
        .def_property_readonly("left_sample",
                               [](const fwdpp::ts::marginal_tree& m) {
                                   return make_1d_ndarray(m.left_sample);
                               })
        .def_property_readonly("right_sample",
                               [](const fwdpp::ts::marginal_tree& m) {
                                   return make_1d_ndarray(m.right_sample);
                               })
        .def_property_readonly("next_sample",
                               [](const fwdpp::ts::marginal_tree& m) {
                                   return make_1d_ndarray(m.next_sample);
                               })
        .def_property_readonly("sample_index_map",
                               [](const fwdpp::ts::marginal_tree& m) {
                                   return make_1d_ndarray(m.sample_index_map);
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
        "Allows left-to-right visiting of the marginal trees embedded in a "
        ":class:`fwdpy11.ts.TableCollection`")
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

    py::class_<VariantIterator>(m, "VariantIterator")
        .def(py::init<const fwdpp::ts::table_collection&,
                      const std::vector<fwdpy11::Mutation>&,
                      const std::vector<fwdpp::ts::TS_NODE_INT>&>(),
             py::keep_alive<1, 2>())
        .def("next_variant", &VariantIterator::next_variant)
        .def("__iter__", [](VariantIterator& v) { return v.next_variant(); },
             py::keep_alive<0, 1>());

    m.def("simplify", &simplify, py::arg("pop"), py::arg("samples"),
          R"delim(
            Simplify a TableCollection.

            :param pop: A :class:`fwdpy11.Population`
            :param samples: A list of samples (node indexes).
                
            :return: The simplified tables and list mapping input sample IDs to output IDS

            :rtype: :class:`fwdpy11.ts.TableCollection`

            Note that the samples argument is agnostic with respect to the time of
            the nodes in the input tables. Thus, you may do things like simplify
            to a set of "currently-alive" nodes plus some or all ancient samples by
            including some node IDs from :attr:`fwdpy11.ts.TableCollection.preserved_nodes`.

            )delim");

    m.def("infinite_sites", [](const fwdpy11::GSLrng_t& rng,
                               fwdpy11::Population& pop, const double mu) {
        std::queue<std::size_t> recycling_bin = fwdpp::ts::make_mut_queue(
            pop.mcounts, pop.mcounts_from_preserved_nodes);
        std::vector<fwdpp::ts::TS_NODE_INT> samples(2 * pop.N);
        std::iota(samples.begin(), samples.end(), 0);
        const auto apply_mutations
            = [&recycling_bin, &rng, &pop, &samples,
               mu](const double left, const double right,
                   const fwdpp::uint_t generation) {
                  return generate_neutral_variants(recycling_bin, pop, rng,
                                                   left, right, generation);
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

    m.def("make_data_matrix",
          [](const fwdpy11::Population& pop,
             const std::vector<fwdpp::ts::TS_NODE_INT>& samples,
             const bool record_neutral, const bool record_selected) {
              return fwdpp::ts::generate_data_matrix(
                  pop.tables, samples, pop.mutations, record_neutral,
                  record_selected);
          },
          py::arg("pop"), py::arg("samples"), py::arg("record_neutral"),
          py::arg("record_selected"),
          R"delim(
     Create a :class:`fwpdy11.sampling.DataMatrix` from a table collection.
     
     :param pop: A population
     :type pop: :class:`fwdpy11.Population`
     :param samples: A list of sample nodes
     :type samples: list
     :param record_neutral: If True, generate data for neutral variants
     :type record_neutral: boolean
     :param record_selected: If True, generate data for selected variants
     :type record_selected: boolean
     )delim");
}
