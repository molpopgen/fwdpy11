#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <fwdpp/ts/std_table_collection.hpp>
#include <fwdpp/ts/tree_visitor.hpp>
#include <fwdpp/ts/marginal_tree.hpp>
#include <fwdpp/ts/marginal_tree_functions/roots.hpp>
#include <fwdpp/ts/marginal_tree_functions/samples.hpp>
#include <fwdpy11/numpy/array.hpp>
#include "node_traversal.hpp"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(fwdpp::ts::std_table_collection::edge_table);
PYBIND11_MAKE_OPAQUE(fwdpp::ts::std_table_collection::node_table);
PYBIND11_MAKE_OPAQUE(fwdpp::ts::std_table_collection::mutation_table);

class tree_visitor_wrapper
{
  private:
    inline fwdpp::ts::table_index_t
    fetch(const std::vector<fwdpp::ts::table_index_t>& data,
          fwdpp::ts::table_index_t index)
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
        if (!std::isfinite(from) || !std::isfinite(until) || from >= genome_length
            || !(until > from))
            {
                throw std::invalid_argument("invalid position range");
            }
    }

    // We hold a reference to the input
    // TableCollection, which prevents it
    // bad things from happening in the
    // calling environment
    std::shared_ptr<fwdpp::ts::std_table_collection> tables_;
    fwdpp::ts::std_table_collection::site_table::const_iterator first_site, end_of_sites,
        current_site;
    fwdpp::ts::std_table_collection::mutation_table::const_iterator first_mutation,
        end_of_mutations, current_mutation;
    bool update_samples;
    const double from, until;

  public:
    fwdpp::ts::tree_visitor<fwdpp::ts::std_table_collection> visitor;
    std::vector<fwdpp::ts::table_index_t> samples_below_buffer;
    tree_visitor_wrapper(std::shared_ptr<fwdpp::ts::std_table_collection> tables,
                         const std::vector<fwdpp::ts::table_index_t>& samples,
                         bool update_samples_below, double start, double stop)
        : tables_(tables), first_site(tables_->sites.begin()),
          end_of_sites(tables_->sites.end()), current_site(first_site),
          first_mutation(tables_->mutations.begin()),
          end_of_mutations(tables_->mutations.end()), current_mutation(first_mutation),
          update_samples(update_samples_below), from(start), until(stop),
          visitor(*tables_, samples,
                  fwdpp::ts::update_samples_list(update_samples_below)),
          samples_below_buffer()
    {
        validate_from_until(tables_->genome_length());
    }

    tree_visitor_wrapper(std::shared_ptr<fwdpp::ts::std_table_collection> tables,
                         const std::vector<fwdpp::ts::table_index_t>& samples,
                         const std::vector<fwdpp::ts::table_index_t>& preserved_nodes,
                         bool update_samples_below, double start, double stop)
        : tables_(tables), first_site(tables_->sites.begin()),
          end_of_sites(tables_->sites.end()), current_site(first_site),
          first_mutation(tables_->mutations.begin()),
          end_of_mutations(tables_->mutations.end()), current_mutation(first_mutation),
          update_samples(update_samples_below), from(start), until(stop),
          visitor(*tables_, samples, preserved_nodes,
                  fwdpp::ts::update_samples_list(update_samples_below)),
          samples_below_buffer()
    {
        validate_from_until(tables_->genome_length());
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
        current_site = std::lower_bound(
            current_site, end_of_sites, pos,
            [](const fwdpp::ts::site& s, double value) { return s.position < value; });
        if (current_site < end_of_sites)
            {
                pos = current_site->position;
                current_mutation = std::lower_bound(
                    current_mutation, end_of_mutations, pos,
                    [this](const fwdpp::ts::mutation_record& mr, double value) {
                        return (first_site + mr.site)->position < value;
                    });
                if (current_mutation < end_of_mutations
                    && (first_site + current_mutation->site)->position != pos)
                    {
                        throw std::runtime_error("error site and mutation iterators");
                    }
            }
        if (visitor.tree().left >= until)
            {
                return false;
            }
        return rv;
    }

    fwdpp::ts::table_index_t
    sample_size() const
    {
        return visitor.tree().sample_size();
    }

    py::array_t<fwdpp::ts::table_index_t>
    nodes()
    {
        std::vector<fwdpp::ts::table_index_t> vnodes(nodes_preorder(visitor.tree()));
        return fwdpy11::make_1d_array_with_capsule(std::move(vnodes));
    }

    py::array_t<fwdpp::ts::table_index_t>
    samples() const
    {
        std::vector<fwdpp::ts::table_index_t> s(visitor.tree().samples_list_begin(),
                                                visitor.tree().samples_list_end());
        return fwdpy11::make_1d_array_with_capsule(std::move(s));
    }

    py::array
    samples_below(const fwdpp::ts::table_index_t node, bool sorted)
    {
        if (!update_samples)
            {
                throw std::invalid_argument("sample tracking not initialized");
            }
        if (node == fwdpp::ts::NULL_INDEX)
            {
                throw std::invalid_argument("invalid node");
            }
        samples_below_buffer.clear();
        fwdpp::ts::process_samples(
            visitor.tree(), fwdpp::ts::convert_sample_index_to_nodes(true), node,
            [this](fwdpp::ts::table_index_t s) { samples_below_buffer.push_back(s); });
        if (sorted)
            {
                std::sort(begin(samples_below_buffer), end(samples_below_buffer));
            }
        return fwdpy11::make_1d_ndarray_readonly(samples_below_buffer);
    }

    fwdpp::ts::table_index_t
    parent(fwdpp::ts::table_index_t u)
    {
        return fetch(this->visitor.tree().parents, u);
    }

    fwdpp::ts::table_index_t
    left_sib(fwdpp::ts::table_index_t u)
    {
        return fetch(this->visitor.tree().left_sib, u);
    }

    fwdpp::ts::table_index_t
    right_sib(fwdpp::ts::table_index_t u)
    {
        return fetch(this->visitor.tree().right_sib, u);
    }

    fwdpp::ts::table_index_t
    left_child(fwdpp::ts::table_index_t u)
    {
        return fetch(this->visitor.tree().left_child, u);
    }

    fwdpp::ts::table_index_t
    right_child(fwdpp::ts::table_index_t u)
    {
        return fetch(this->visitor.tree().right_child, u);
    }

    fwdpp::ts::table_index_t
    leaf_counts(fwdpp::ts::table_index_t u)
    {
        return fetch(this->visitor.tree().leaf_counts, u);
    }

    fwdpp::ts::table_index_t
    preserved_leaf_counts(fwdpp::ts::table_index_t u)
    {
        return fetch(this->visitor.tree().preserved_leaf_counts, u);
    }

    std::shared_ptr<fwdpp::ts::std_table_collection>
    get_tables() const
    {
        return tables_;
    }

    std::pair<fwdpp::ts::std_table_collection::site_table::const_iterator,
              fwdpp::ts::std_table_collection::site_table::const_iterator>
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
        auto end_of_range = std::lower_bound(
            current_site, end_of_sites, pos,
            [](const fwdpp::ts::site& s, double value) { return s.position < value; });
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

    std::pair<fwdpp::ts::std_table_collection::mutation_table::const_iterator,
              fwdpp::ts::std_table_collection::mutation_table::const_iterator>
    get_mutations_on_current_tree()
    {
        double pos = std::min(visitor.tree().right, until);
        while ((current_site < end_of_sites)
               && (current_site->position < visitor.tree().left))
            {
                ++current_site;
            }
        while ((current_site < end_of_sites) && (current_mutation < end_of_mutations)
               && (first_site + current_mutation->site)->position
                      != current_site->position)
            {
                ++current_mutation;
            }
        if ((current_site == end_of_sites) || (current_mutation == end_of_mutations)
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
    py::class_<tree_visitor_wrapper>(m, "ll_TreeIterator")
        .def(py::init<std::shared_ptr<fwdpp::ts::std_table_collection>,
                      const std::vector<fwdpp::ts::table_index_t>&, bool, double,
                      double>(),
             py::arg("tables"), py::arg("samples"), py::arg("update_samples") = false,
             py::arg("begin") = 0.0, py::arg("end") = std::numeric_limits<double>::max())
        .def(py::init<std::shared_ptr<fwdpp::ts::std_table_collection>,
                      const std::vector<fwdpp::ts::table_index_t>&,
                      const std::vector<fwdpp::ts::table_index_t>&, bool, double,
                      double>(),
             py::arg("tables"), py::arg("samples"), py::arg("ancient_samples"),
             py::arg("update_samples") = false, py::arg("begin") = 0.0,
             py::arg("end") = std::numeric_limits<double>::max())
        .def("_parent", &tree_visitor_wrapper::parent, "Return parent of a node")
        .def("_leaf_counts", &tree_visitor_wrapper::leaf_counts)
        .def("_preserved_leaf_counts", &tree_visitor_wrapper::preserved_leaf_counts)
        .def("_left_sib", &tree_visitor_wrapper::left_sib)
        .def("_right_sib", &tree_visitor_wrapper::right_sib)
        .def("_left_child", &tree_visitor_wrapper::left_child)
        .def("_right_child", &tree_visitor_wrapper::right_child)
        .def_property_readonly(
            "_left",
            [](const tree_visitor_wrapper& self) { return self.visitor.tree().left; })
        .def_property_readonly(
            "_right",
            [](const tree_visitor_wrapper& self) { return self.visitor.tree().right; })
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
             [](tree_visitor_wrapper& self) -> tree_visitor_wrapper& { return self; })
        .def("_total_time",
             [](const tree_visitor_wrapper& self,
                const fwdpp::ts::std_table_collection::node_table& nodes) {
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
                         if (m.parents[i] != fwdpp::ts::NULL_INDEX)
                             {
                                 tt += nodes[i].time - nodes[m.parents[i]].time;
                             }
                     }
                 return tt;
             })
        .def_property_readonly("sample_size", &tree_visitor_wrapper::sample_size)
        .def_property_readonly(
            "_roots",
            [](const tree_visitor_wrapper& self) {
                auto roots = fwdpp::ts::get_roots(self.visitor.tree());
                return fwdpy11::make_1d_array_with_capsule(std::move(roots));
            })
        .def("_nodes", &tree_visitor_wrapper::nodes)
        .def("_samples", &tree_visitor_wrapper::samples)
        .def_property_readonly("_tables", &tree_visitor_wrapper::get_tables)
        .def("_samples_below", &tree_visitor_wrapper::samples_below, py::arg("node"),
             py::arg("sorted") = false)
        .def(
            "_sites",
            [](tree_visitor_wrapper& self) {
                auto rv = self.get_sites_on_current_tree();
                return py::make_iterator(rv.first, rv.second);
            },
            py::keep_alive<0, 1>())
        .def(
            "_mutations",
            [](tree_visitor_wrapper& self) {
                auto r = self.get_mutations_on_current_tree();
                return py::make_iterator(r.first, r.second);
            },
            py::keep_alive<0, 1>());
}
