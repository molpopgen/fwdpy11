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

    bool update_samples;
    const double from, until;

  public:
    fwdpp::ts::tree_visitor visitor;
    std::vector<fwdpp::ts::TS_NODE_INT> samples_below_buffer;
    tree_visitor_wrapper(const fwdpp::ts::table_collection& tables,
                         const std::vector<fwdpp::ts::TS_NODE_INT>& samples,
                         bool update_samples_below, double start, double stop)
        : update_samples(update_samples_below), from(start), until(stop),
          visitor(tables, samples,
                  fwdpp::ts::update_samples_list(update_samples_below)),
          samples_below_buffer()
    {
        validate_from_until(tables.genome_length());
    }

    tree_visitor_wrapper(
        const fwdpp::ts::table_collection& tables,
        const std::vector<fwdpp::ts::TS_NODE_INT>& samples,
        const std::vector<fwdpp::ts::TS_NODE_INT>& preserved_nodes,
        bool update_samples_below, double start, double stop)
        : update_samples(update_samples_below), from(start), until(stop),
          visitor(tables, samples,
                  fwdpp::ts::update_samples_list(update_samples_below)),
          samples_below_buffer()
    {
        validate_from_until(tables.genome_length());
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
        std::vector<fwdpp::ts::TS_NODE_INT>* nodes
            = new std::vector<fwdpp::ts::TS_NODE_INT>(
                nodes_preorder(visitor.tree()));

        auto capsule = py::capsule(nodes, [](void* x) {
            delete reinterpret_cast<std::vector<fwdpp::ts::TS_NODE_INT>*>(x);
        });
        return py::array_t<fwdpp::ts::TS_NODE_INT>(nodes->size(),
                                                   nodes->data(), capsule);
    }

    py::array
    samples() const
    {
        std::vector<fwdpp::ts::TS_NODE_INT>* s
            = new std::vector<fwdpp::ts::TS_NODE_INT>(
                visitor.tree().samples_list_begin(),
                visitor.tree().samples_list_end());
        auto capsule = py::capsule(s, [](void* x) {
            delete reinterpret_cast<std::vector<fwdpp::ts::TS_NODE_INT>*>(x);
        });
        return py::array(s->size(), s->data(), capsule);
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
        auto capsule = py::capsule(&samples_below_buffer, [](void* x) {
            reinterpret_cast<decltype(samples_below_buffer)*>(x)->clear();
        });
        return py::array(samples_below_buffer.size(),
                         samples_below_buffer.data(), capsule);
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
        .def(py::init<const fwdpp::ts::table_collection&,
                      const std::vector<fwdpp::ts::TS_NODE_INT>&, bool, double,
                      double>(),
             py::arg("tables"), py::arg("samples"),
             py::arg("update_samples") = false, py::arg("begin") = 0.0,
             py::arg("end") = std::numeric_limits<double>::max())
        .def(py::init<const fwdpp::ts::table_collection&,
                      const std::vector<fwdpp::ts::TS_NODE_INT>&,
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
                std::vector<fwdpp::ts::TS_NODE_INT>* roots
                    = new std::vector<fwdpp::ts::TS_NODE_INT>(
                        fwdpp::ts::get_roots(self.visitor.tree()));

                auto capsule = py::capsule(roots, [](void* x) {
                    delete reinterpret_cast<
                        std::vector<fwdpp::ts::TS_NODE_INT>*>(x);
                });

                return py::array(roots->size(), roots->data(), capsule);
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
             py::arg("node"), py::arg("sorted") = false);
}
