#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <fwdpp/ts/tree_visitor.hpp>
#include <fwdpp/ts/marginal_tree.hpp>

namespace py = pybind11;

template <typename T>
inline py::array_t<T>
make_1d_ndarray(const std::vector<T>& v)
{
    auto rv
        = py::array_t<T>({ v.size() }, { sizeof(T) }, v.data(), py::cast(v));
    return rv;
}

PYBIND11_MAKE_OPAQUE(fwdpp::ts::edge_vector);
PYBIND11_MAKE_OPAQUE(fwdpp::ts::node_vector);
PYBIND11_MAKE_OPAQUE(fwdpp::ts::mutation_key_vector);


struct tree_visitor_wrapper
{
    fwdpp::ts::tree_visitor visitor;
    py::array parents, leaf_counts, preserved_leaf_counts, left_sib, right_sib,
        left_child, right_child, left_sample, right_sample, next_sample,
        sample_index_map;
    tree_visitor_wrapper(const fwdpp::ts::table_collection& tables,
                         const std::vector<fwdpp::ts::TS_NODE_INT>& samples)
        : visitor(tables, samples),
          parents(make_1d_ndarray(visitor.tree().parents)),
          leaf_counts(make_1d_ndarray(visitor.tree().leaf_counts)),
          preserved_leaf_counts(
              make_1d_ndarray(visitor.tree().preserved_leaf_counts)),
          left_sib(make_1d_ndarray(visitor.tree().left_sib)),
          right_sib(make_1d_ndarray(visitor.tree().right_sib)),
          left_child(make_1d_ndarray(visitor.tree().left_child)),
          right_child(make_1d_ndarray(visitor.tree().right_child)),
          left_sample(make_1d_ndarray(visitor.tree().left_sample)),
          right_sample(make_1d_ndarray(visitor.tree().right_sample)),
          next_sample(make_1d_ndarray(visitor.tree().next_sample)),
          sample_index_map(make_1d_ndarray(visitor.tree().sample_index_map))
    {
    }

    inline bool
    operator()(const bool update_samples)
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
                      const std::vector<fwdpp::ts::TS_NODE_INT>&>(),
             py::arg("tables"), py::arg("samples"))
        .def_readwrite("parents", &tree_visitor_wrapper::parents)
        .def_readonly("leaf_counts", &tree_visitor_wrapper::leaf_counts)
        .def_readonly("preserved_leaf_counts",
                      &tree_visitor_wrapper::preserved_leaf_counts)
        .def_readonly("left_sib", &tree_visitor_wrapper::left_sib)
        .def_readonly("right_sib", &tree_visitor_wrapper::right_sib)
        .def_readonly("left_child", &tree_visitor_wrapper::left_child)
        .def_readonly("right_child", &tree_visitor_wrapper::right_child)
        .def_readonly("left_sample", &tree_visitor_wrapper::leaf_counts)
        .def_readonly("right_sample", &tree_visitor_wrapper::leaf_counts)
        .def_readonly("next_sample", &tree_visitor_wrapper::next_sample)
        .def_readonly("sample_index_map",
                      &tree_visitor_wrapper::sample_index_map)
        .def_property_readonly("left",
                               [](const tree_visitor_wrapper& self) {
                                   return self.visitor.tree().left;
                               })
        .def_property_readonly("right",
                               [](const tree_visitor_wrapper& self) {
                                   return self.visitor.tree().right;
                               })
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
                               [](const tree_visitor_wrapper& self) {
                                   return self.sample_size();
                               });
}
