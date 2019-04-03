#include <fwdpp/ts/mutation_record.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <sstream>

namespace py = pybind11;

void
init_ts_MutationRecord(py::module& m)
{
    PYBIND11_NUMPY_DTYPE(fwdpp::ts::mutation_record, node, key);
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
                      "Index of the mutation in the population")
        .def("__repr__",
             [](const fwdpp::ts::mutation_record& self) {
                 std::ostringstream out;
                 out << "MutationRecord(node=" << self.node
                     << ", key=" << self.key << ")";
                 return out.str();
             })
        .def(py::pickle(
            [](const fwdpp::ts::mutation_record& m) {
                return py::make_tuple(m.node, m.key);
            },
            [](py::tuple t) {
                return fwdpp::ts::mutation_record{
                    t[0].cast<decltype(fwdpp::ts::mutation_record::node)>(),
                    t[1].cast<decltype(fwdpp::ts::mutation_record::key)>()
                };
            }));
}

