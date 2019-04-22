#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <fwdpp/ts/tree_visitor.hpp>
#include <fwdpp/ts/marginal_tree.hpp>
#include <fwdpy11/numpy/array.hpp>

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
    std::vector<fwdpp::ts::TS_NODE_INT> sample_list_buffer;
    tree_visitor_wrapper(const fwdpp::ts::table_collection& tables,
                         const std::vector<fwdpp::ts::TS_NODE_INT>& samples,
                         bool update_sample_list, double start, double stop)
        : update_samples(update_sample_list), from(start), until(stop),
          visitor(tables, samples), sample_list_buffer()
    {
        validate_from_until(tables.genome_length());
    }

    tree_visitor_wrapper(
        const fwdpp::ts::table_collection& tables,
        const std::vector<fwdpp::ts::TS_NODE_INT>& samples,
        const std::vector<fwdpp::ts::TS_NODE_INT>& preserved_nodes,
        bool update_sample_list, double start, double stop)
        : update_samples(update_sample_list), from(start), until(stop),
          visitor(tables, samples), sample_list_buffer()
    {
        validate_from_until(tables.genome_length());
    }

    inline bool
    operator()()
    {
        bool rv = false;
        if (update_samples)
            {
                rv = visitor(std::true_type(), std::true_type());
                while (visitor.tree().right < from)
                    {
                        rv = visitor(std::true_type(), std::true_type());
                    }
            }
        else
            {
                rv = visitor(std::true_type(), std::false_type());
                while (visitor.tree().right < from)
                    {
                        rv = visitor(std::true_type(), std::false_type());
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
        return visitor.tree().sample_size;
    }

    py::array
    sample_list(const fwdpp::ts::TS_NODE_INT node, bool sorted)
    {
        if (!update_samples)
            {
                throw std::invalid_argument("sample tracking not initialized");
            }
        if (node == fwdpp::ts::TS_NULL_NODE)
            {
                throw std::invalid_argument("invalid node");
            }
        sample_list_buffer.clear();
        const auto& marginal = visitor.tree();
        auto right = marginal.right_sample[node];
        auto index = marginal.left_sample[node];
        while (true)
            {
                sample_list_buffer.push_back(index);
                if (index == right)
                    {
                        break;
                    }
                index = marginal.next_sample[index];
            }
        if (sorted)
            {
                std::sort(begin(sample_list_buffer), end(sample_list_buffer));
            }
        auto capsule = py::capsule(&sample_list_buffer, [](void* x) {
            reinterpret_cast<decltype(sample_list_buffer)*>(x)->clear();
        });
        return py::array(sample_list_buffer.size(), sample_list_buffer.data(),
                         capsule);
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
             py::arg("update_sample_list") = false, py::arg("begin") = 0.0,
             py::arg("end") = std::numeric_limits<double>::max())
        .def(py::init<const fwdpp::ts::table_collection&,
                      const std::vector<fwdpp::ts::TS_NODE_INT>&,
                      const std::vector<fwdpp::ts::TS_NODE_INT>&, bool, double,
                      double>(),
             py::arg("tables"), py::arg("samples"), py::arg("ancient_samples"),
             py::arg("update_sample_list") = false, py::arg("begin") = 0.0,
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
                    = new std::vector<fwdpp::ts::TS_NODE_INT>();

                auto capsule = py::capsule(roots, [](void* x) {
                    delete reinterpret_cast<
                        std::vector<fwdpp::ts::TS_NODE_INT>*>(x);
                });

                return py::array(roots->size(), roots->data(), capsule);
            },
            R"delim(
            Return marginal tree roots as numpy.ndarray
            
            .. versionadded:: 0.4.0
            )delim")
        .def("sample_list", &tree_visitor_wrapper::sample_list,
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
