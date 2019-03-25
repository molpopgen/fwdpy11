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

struct tree_visitor_wrapper
{
    fwdpp::ts::tree_visitor visitor;
    py::array parents, leaf_counts, preserved_leaf_counts, left_sib, right_sib,
        left_child, right_child, left_sample, right_sample, next_sample,
        sample_index_map;
    bool update_samples;
    std::vector<fwdpp::ts::TS_NODE_INT> sample_list_buffer;
    tree_visitor_wrapper(const fwdpp::ts::table_collection& tables,
                         const std::vector<fwdpp::ts::TS_NODE_INT>& samples,
                         bool update_sample_list)
        : visitor(tables, samples),
          parents(fwdpy11::make_1d_ndarray(visitor.tree().parents)),
          leaf_counts(fwdpy11::make_1d_ndarray(visitor.tree().leaf_counts)),
          preserved_leaf_counts(
              fwdpy11::make_1d_ndarray(visitor.tree().preserved_leaf_counts)),
          left_sib(fwdpy11::make_1d_ndarray(visitor.tree().left_sib)),
          right_sib(fwdpy11::make_1d_ndarray(visitor.tree().right_sib)),
          left_child(fwdpy11::make_1d_ndarray(visitor.tree().left_child)),
          right_child(fwdpy11::make_1d_ndarray(visitor.tree().right_child)),
          left_sample(fwdpy11::make_1d_ndarray(visitor.tree().left_sample)),
          right_sample(fwdpy11::make_1d_ndarray(visitor.tree().right_sample)),
          next_sample(fwdpy11::make_1d_ndarray(visitor.tree().next_sample)),
          sample_index_map(
              fwdpy11::make_1d_ndarray(visitor.tree().sample_index_map)),
          update_samples(update_sample_list), sample_list_buffer()
    {
    }

    tree_visitor_wrapper(
        const fwdpp::ts::table_collection& tables,
        const std::vector<fwdpp::ts::TS_NODE_INT>& samples,
        const std::vector<fwdpp::ts::TS_NODE_INT>& preserved_nodes,
        bool update_sample_list)
        : visitor(tables, samples),
          parents(fwdpy11::make_1d_ndarray(visitor.tree().parents)),
          leaf_counts(fwdpy11::make_1d_ndarray(visitor.tree().leaf_counts)),
          preserved_leaf_counts(
              fwdpy11::make_1d_ndarray(visitor.tree().preserved_leaf_counts)),
          left_sib(fwdpy11::make_1d_ndarray(visitor.tree().left_sib)),
          right_sib(fwdpy11::make_1d_ndarray(visitor.tree().right_sib)),
          left_child(fwdpy11::make_1d_ndarray(visitor.tree().left_child)),
          right_child(fwdpy11::make_1d_ndarray(visitor.tree().right_child)),
          left_sample(fwdpy11::make_1d_ndarray(visitor.tree().left_sample)),
          right_sample(fwdpy11::make_1d_ndarray(visitor.tree().right_sample)),
          next_sample(fwdpy11::make_1d_ndarray(visitor.tree().next_sample)),
          sample_index_map(
              fwdpy11::make_1d_ndarray(visitor.tree().sample_index_map)),
          update_samples(update_sample_list), sample_list_buffer()
    {
    }

    inline bool
    operator()()
    {
        bool rv = false;
        if (update_samples)
            {
                rv = visitor(std::true_type(), std::true_type());
            }
        else
            {
                rv = visitor(std::true_type(), std::false_type());
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
};

void
init_tree_iterator(py::module& m)
{
    py::class_<tree_visitor_wrapper>(m, "TreeIterator",
                                     R"delim(
            Iterate over the marginal trees in a :class:`fwdpy11.ts.TableCollection`
            
            .. versionadded 0.3.0
            )delim")
        .def(py::init<const fwdpp::ts::table_collection&,
                      const std::vector<fwdpp::ts::TS_NODE_INT>&, bool>(),
             py::arg("tables"), py::arg("samples"),
             py::arg("update_sample_list") = false)
        .def(py::init<const fwdpp::ts::table_collection&,
                      const std::vector<fwdpp::ts::TS_NODE_INT>&,
                      const std::vector<fwdpp::ts::TS_NODE_INT>&, bool>(),
             py::arg("tables"), py::arg("samples"), py::arg("ancient_samples"),
             py::arg("update_sample_list") = false)
        .def_readwrite("parents", &tree_visitor_wrapper::parents,
                       "Vector of child -> parent relationships")
        .def_readonly("leaf_counts", &tree_visitor_wrapper::leaf_counts,
                      "Leaf counts for each node")
        .def_readonly("preserved_leaf_counts",
                      &tree_visitor_wrapper::preserved_leaf_counts,
                      "Ancient sample leaf counts for each node")
        .def_readonly("left_sib", &tree_visitor_wrapper::left_sib,
                      "Return the left sibling of the current node")
        .def_readonly("right_sib", &tree_visitor_wrapper::right_sib,
                      "Return the right sibling of the current node")
        .def_readonly("left_child", &tree_visitor_wrapper::left_child,
                      "Mapping of current node id to its left child")
        .def_readonly("right_child", &tree_visitor_wrapper::right_child,
                      "Mapping of current node id to its right child")
        .def_readonly("left_sample", &tree_visitor_wrapper::left_sample)
        .def_readonly("right_sample", &tree_visitor_wrapper::right_sample)
        .def_readonly("next_sample", &tree_visitor_wrapper::next_sample)
        .def_readonly("sample_index_map",
                      &tree_visitor_wrapper::sample_index_map)
        .def_property_readonly("left",
                               [](const tree_visitor_wrapper& self) {
                                   return self.visitor.tree().left;
                               },
                               "Left edge of genomic interval (inclusive)")
        .def_property_readonly("right",
                               [](const tree_visitor_wrapper& self) {
                                   return self.visitor.tree().right;
                               },
                               "Right edge of genomic interval (exclusive)")
        .def("__call__", &tree_visitor_wrapper::operator())
        .def("total_time",
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
                                 tt += nodes[i].time
                                       - nodes[m.parents[i]].time;
                             }
                     }
                 return tt;
             },
             "Return the sum of branch lengths")
        .def_property_readonly("sample_size",
                               &tree_visitor_wrapper::sample_size)
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
