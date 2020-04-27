#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <fwdpp/ts/tree_visitor.hpp>
#include <fwdpp/ts/marginal_tree.hpp>
#include <fwdpp/ts/marginal_tree_functions/roots.hpp>
#include <fwdpp/ts/marginal_tree_functions/samples.hpp>
#include <fwdpy11/numpy/array.hpp>
#include "node_traversal.hpp"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(fwdpp::ts::edge_vector);
PYBIND11_MAKE_OPAQUE(fwdpp::ts::node_vector);
PYBIND11_MAKE_OPAQUE(fwdpp::ts::mutation_key_vector);

class tree_visitor_wrapper
{
  private:
    inline fwdpp::ts::TS_NODE_INT
    fetch(const std::vector<fwdpp::ts::TS_NODE_INT>& data,
          fwdpp::ts::TS_NODE_INT index)
    {
        if (static_cast<std::size_t>(index) >= data.size())
            {
                throw std::out_of_range("node index is our of range");
            }
        return data[index];
    }

    void
    validate_from_until(const double genome_length)
    {
        if (!std::isfinite(from) || !std::isfinite(until)
            || from >= genome_length || !(until > from))
            {
                throw std::invalid_argument("invalid position range");
            }
    }

    // We hold a reference to the input
    // TableCollection, which prevents it
    // bad things from happening in the
    // calling environment
    py::object tables_;
    fwdpp::ts::site_vector::const_iterator first_site, end_of_sites,
        current_site;
    fwdpp::ts::mutation_key_vector::const_iterator first_mutation,
        end_of_mutations, current_mutation;
    bool update_samples;
    const double from, until;

  public:
    fwdpp::ts::tree_visitor visitor;
    std::vector<fwdpp::ts::TS_NODE_INT> samples_below_buffer;
    tree_visitor_wrapper(py::object tables,
                         const std::vector<fwdpp::ts::TS_NODE_INT>& samples,
                         bool update_samples_below, double start, double stop)
        : tables_(tables),
          first_site(tables_.cast<const fwdpp::ts::table_collection&>()
                         .site_table.begin()),
          end_of_sites(tables_.cast<const fwdpp::ts::table_collection&>()
                           .site_table.end()),
          current_site(first_site),
          first_mutation(tables_.cast<const fwdpp::ts::table_collection&>()
                             .mutation_table.begin()),
          end_of_mutations(tables_.cast<const fwdpp::ts::table_collection&>()
                               .mutation_table.end()),
          current_mutation(first_mutation),
          update_samples(update_samples_below), from(start), until(stop),
          visitor(tables_.cast<const fwdpp::ts::table_collection&>(), samples,
                  fwdpp::ts::update_samples_list(update_samples_below)),
          samples_below_buffer()
    {
        validate_from_until(
            tables_.cast<fwdpp::ts::table_collection&>().genome_length());
    }

    tree_visitor_wrapper(
        py::object tables, const std::vector<fwdpp::ts::TS_NODE_INT>& samples,
        const std::vector<fwdpp::ts::TS_NODE_INT>& preserved_nodes,
        bool update_samples_below, double start, double stop)
        : tables_(tables),
          first_site(tables_.cast<const fwdpp::ts::table_collection&>()
                         .site_table.begin()),
          end_of_sites(tables_.cast<const fwdpp::ts::table_collection&>()
                           .site_table.end()),
          current_site(first_site),
          first_mutation(tables_.cast<const fwdpp::ts::table_collection&>()
                             .mutation_table.begin()),
          end_of_mutations(tables_.cast<const fwdpp::ts::table_collection&>()
                               .mutation_table.end()),
          current_mutation(first_mutation),
          update_samples(update_samples_below), from(start), until(stop),
          visitor(tables_.cast<const fwdpp::ts::table_collection&>(), samples,
                  fwdpp::ts::update_samples_list(update_samples_below)),
          samples_below_buffer()
    {
        validate_from_until(
            tables_.cast<fwdpp::ts::table_collection&>().genome_length());
    }

    inline bool
    operator()()
    {
        bool rv = false;
        rv = visitor();
        while (visitor.tree().right < from)
            {
                rv = visitor();
            }
        double pos = std::max(visitor.tree().left, from);
        current_site
            = std::lower_bound(current_site, end_of_sites, pos,
                               [](const fwdpp::ts::site& s, double value) {
                                   return s.position < value;
                               });
        if (current_site < end_of_sites)
            {
                pos = current_site->position;
                current_mutation = std::lower_bound(
                    current_mutation, end_of_mutations, pos,
                    [this](const fwdpp::ts::mutation_record& mr,
                           double value) {
                        return (first_site + mr.site)->position < value;
                    });
                if (current_mutation < end_of_mutations
                    && (first_site + current_mutation->site)->position != pos)
                    {
                        throw std::runtime_error(
                            "error site and mutation iterators");
                    }
            }
        if (visitor.tree().left >= until)
            {
                return false;
            }
        return rv;
    }

    fwdpp::ts::TS_NODE_INT
    sample_size() const
    {
        return visitor.tree().sample_size();
    }

    py::array_t<fwdpp::ts::TS_NODE_INT>
    nodes()
    {
        std::vector<fwdpp::ts::TS_NODE_INT> vnodes(
            nodes_preorder(visitor.tree()));
        return fwdpy11::make_1d_array_with_capsule(std::move(vnodes));
    }

    py::array_t<fwdpp::ts::TS_NODE_INT>
    samples() const
    {
        std::vector<fwdpp::ts::TS_NODE_INT> s(
            visitor.tree().samples_list_begin(),
            visitor.tree().samples_list_end());
        return fwdpy11::make_1d_array_with_capsule(std::move(s));
    }

    py::array
    samples_below(const fwdpp::ts::TS_NODE_INT node, bool sorted)
    {
        if (!update_samples)
            {
                throw std::invalid_argument("sample tracking not initialized");
            }
        if (node == fwdpp::ts::TS_NULL_NODE)
            {
                throw std::invalid_argument("invalid node");
            }
        samples_below_buffer.clear();
        fwdpp::ts::process_samples(
            visitor.tree(), fwdpp::ts::convert_sample_index_to_nodes(true),
            node, [this](fwdpp::ts::TS_NODE_INT s) {
                samples_below_buffer.push_back(s);
            });
        if (sorted)
            {
                std::sort(begin(samples_below_buffer),
                          end(samples_below_buffer));
            }
        return fwdpy11::make_1d_ndarray_readonly(samples_below_buffer);
    }

    fwdpp::ts::TS_NODE_INT
    parent(fwdpp::ts::TS_NODE_INT u)
    {
        return fetch(this->visitor.tree().parents, u);
    }

    fwdpp::ts::TS_NODE_INT
    left_sib(fwdpp::ts::TS_NODE_INT u)
    {
        return fetch(this->visitor.tree().left_sib, u);
    }

    fwdpp::ts::TS_NODE_INT
    right_sib(fwdpp::ts::TS_NODE_INT u)
    {
        return fetch(this->visitor.tree().right_sib, u);
    }

    fwdpp::ts::TS_NODE_INT
    left_child(fwdpp::ts::TS_NODE_INT u)
    {
        return fetch(this->visitor.tree().left_child, u);
    }

    fwdpp::ts::TS_NODE_INT
    right_child(fwdpp::ts::TS_NODE_INT u)
    {
        return fetch(this->visitor.tree().right_child, u);
    }

    fwdpp::ts::TS_NODE_INT
    leaf_counts(fwdpp::ts::TS_NODE_INT u)
    {
        return fetch(this->visitor.tree().leaf_counts, u);
    }

    fwdpp::ts::TS_NODE_INT
    preserved_leaf_counts(fwdpp::ts::TS_NODE_INT u)
    {
        return fetch(this->visitor.tree().preserved_leaf_counts, u);
    }

    py::object
    get_tables() const
    {
        return tables_;
    }

    std::pair<fwdpp::ts::site_vector::const_iterator,
              fwdpp::ts::site_vector::const_iterator>
    get_sites_on_current_tree()
    {
        double pos = std::min(visitor.tree().right, until);
        while (current_site->position < visitor.tree().left)
            {
                ++current_site;
            }
        if ((current_site == end_of_sites)
            || (current_site < end_of_sites && current_site->position >= pos))
            {
                // Return an empty range
                return std::make_pair(end_of_sites, end_of_sites);
            }
        // Find first Site > current_tree.right
        auto end_of_range
            = std::lower_bound(current_site, end_of_sites, pos,
                               [](const fwdpp::ts::site& s, double value) {
                                   return s.position < value;
                               });
        // NOTE: this guards against a corner case.  If the first mutation
        // not in this tree has position == pos, then upper_bound would result
        // in it being included in he interval, which is bad. So, we instead
        // search using lower_bound and check:
        if (end_of_range < end_of_sites && end_of_range->position == pos)
            {
                ++end_of_range;
            }
        return std::make_pair(current_site, end_of_range);
    }

    std::pair<fwdpp::ts::mutation_key_vector::const_iterator,
              fwdpp::ts::mutation_key_vector::const_iterator>
    get_mutations_on_current_tree()
    {
        double pos = std::min(visitor.tree().right, until);
        while ((current_site < end_of_sites)
               && (current_site->position < visitor.tree().left))
            {
                ++current_site;
            }
        while ((current_site < end_of_sites)
               && (current_mutation < end_of_mutations)
               && (first_site + current_mutation->site)->position
                      != current_site->position)
            {
                ++current_mutation;
            }
        if ((current_site == end_of_sites)
            || (current_mutation == end_of_mutations)
            || (current_site < end_of_sites && current_site->position >= pos))
            {
                // Return an empty range
                return std::make_pair(end_of_mutations, end_of_mutations);
            }
        auto end_of_range = std::lower_bound(
            current_mutation, end_of_mutations, pos,
            [this](const fwdpp::ts::mutation_record& mr, double value) {
                return (first_site + mr.site)->position < value;
            });
        if (end_of_range < end_of_mutations
            && (first_site + end_of_range->site)->position == pos)
            {
                ++end_of_range;
            }
        return std::make_pair(current_mutation, end_of_range);
    }
};

void
init_tree_iterator(py::module& m)
{
    py::class_<tree_visitor_wrapper>(m, "TreeIterator",
                                     R"delim(
            Iterate over the marginal trees in a :class:`fwdpy11.TableCollection`
            
            .. versionadded 0.3.0

            .. versionchanged:: 0.4.1
        
                Add begin, end options as floats for initializing
            )delim")
        .def(py::init<py::object, const std::vector<fwdpp::ts::TS_NODE_INT>&,
                      bool, double, double>(),
             py::arg("tables"), py::arg("samples"),
             py::arg("update_samples") = false, py::arg("begin") = 0.0,
             py::arg("end") = std::numeric_limits<double>::max())
        .def(py::init<py::object, const std::vector<fwdpp::ts::TS_NODE_INT>&,
                      const std::vector<fwdpp::ts::TS_NODE_INT>&, bool, double,
                      double>(),
             py::arg("tables"), py::arg("samples"), py::arg("ancient_samples"),
             py::arg("update_samples") = false, py::arg("begin") = 0.0,
             py::arg("end") = std::numeric_limits<double>::max())
        .def("parent", &tree_visitor_wrapper::parent,
             "Return parent of a node")
        .def("leaf_counts", &tree_visitor_wrapper::leaf_counts,
             "Leaf counts for a node")
        .def("preserved_leaf_counts",
             &tree_visitor_wrapper::preserved_leaf_counts,
             "Ancient sample leaf counts for a node")
        .def("left_sib", &tree_visitor_wrapper::left_sib,
             "Return the left sibling of the current node")
        .def("right_sib", &tree_visitor_wrapper::right_sib,
             "Return the right sibling of the current node")
        .def("left_child", &tree_visitor_wrapper::left_child,
             "Mapping of current node id to its left child")
        .def("right_child", &tree_visitor_wrapper::right_child,
             "Mapping of current node id to its right child")
        .def_property_readonly(
            "left",
            [](const tree_visitor_wrapper& self) {
                return self.visitor.tree().left;
            },
            "Left edge of genomic interval (inclusive)")
        .def_property_readonly(
            "right",
            [](const tree_visitor_wrapper& self) {
                return self.visitor.tree().right;
            },
            "Right edge of genomic interval (exclusive)")
        .def("__next__",
             [](tree_visitor_wrapper& self) -> tree_visitor_wrapper& {
                 auto x = self();
                 if (x == false)
                     {
                         throw py::stop_iteration();
                     }
                 return self;
             })
        .def("__iter__",
             [](tree_visitor_wrapper& self) -> tree_visitor_wrapper& {
                 return self;
             })
        .def(
            "total_time",
            [](const tree_visitor_wrapper& self,
               const fwdpp::ts::node_vector& nodes) {
                const auto& m = self.visitor.tree();
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
                                tt += nodes[i].time - nodes[m.parents[i]].time;
                            }
                    }
                return tt;
            },
            "Return the sum of branch lengths")
        .def_property_readonly("sample_size",
                               &tree_visitor_wrapper::sample_size)
        .def_property_readonly(
            "roots",
            [](const tree_visitor_wrapper& self) {
                auto roots = fwdpp::ts::get_roots(self.visitor.tree());
                return fwdpy11::make_1d_array_with_capsule(std::move(roots));
            },
            R"delim(
            Return marginal tree roots as numpy.ndarray
            
            .. versionadded:: 0.4.0

            .. versionchanged:: 0.5.1

                Fixed to not return an empty array.
            )delim")
        .def("nodes", &tree_visitor_wrapper::nodes,
             R"delim("Return the nodes in the current tree.
             
             The return order is preorder.
             
             .. versionadded:: 0.5.1
             )delim")
        .def("samples", &tree_visitor_wrapper::samples,
             "Return the complete sample list")
        .def_property_readonly("tables", &tree_visitor_wrapper::get_tables,
                               "Return the TableCollection")
        .def("samples_below", &tree_visitor_wrapper::samples_below,
             R"delim(
            Return the list of samples descending from a node.

            :param node: A node id
            :type node: int
            :param sorted: (False) Whether or not to sort sample node IDs.
            :type sorted: boolean

            .. note::

                Do not store these sample lists without making a "deep"
                copy.  The internal buffer is re-used.
            )delim",
             py::arg("node"), py::arg("sorted") = false)
        .def(
            "sites",
            [](tree_visitor_wrapper& self) {
                auto rv = self.get_sites_on_current_tree();
                return py::make_iterator(rv.first, rv.second);
            },
            py::keep_alive<0, 1>(),
            R"delim(
            Return iterator over all :class:`fwdpy11.Site` objects
            on the current tree.

            .. versionadded:: 0.5.1
            )delim")
        .def(
            "mutations",
            [](tree_visitor_wrapper& self) {
                auto r = self.get_mutations_on_current_tree();
                return py::make_iterator(r.first, r.second);
            },
            py::keep_alive<0, 1>(),
            R"delim(
            Return iterator over all :class:`fwdpy11.MutationRecord` objects
            on the current tree.

            .. versionadded:: 0.5.1
            )delim");
}
