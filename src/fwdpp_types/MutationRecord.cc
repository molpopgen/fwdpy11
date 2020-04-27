#include <fwdpp/ts/mutation_record.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <sstream>

namespace py = pybind11;

void
init_ts_MutationRecord(py::module& m)
{
    PYBIND11_NUMPY_DTYPE(fwdpp::ts::mutation_record, node, key, site,
                         derived_state, neutral);
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
        .def_readonly("site", &fwdpp::ts::mutation_record::site,
                      R"delim(Index of the mutation's site in the population.
                      
                      .. versionadded:: 0.5.0)delim")
        .def_readonly("derived_state",
                      &fwdpp::ts::mutation_record::derived_state,
                      R"delim(The derived state of the mutation.
                      
                      .. versionadded:: 0.5.0)delim")
        .def_readonly("neutral", &fwdpp::ts::mutation_record::neutral,
                      R"delim(Boolean descriptor of whether or not the mutation
                      affects fitness.
                      
                      .. versionadded:: 0.5.0)delim")
        .def("__repr__",
             [](const fwdpp::ts::mutation_record& self) {
                 std::ostringstream out;
                 out << "MutationRecord(node=" << self.node
                     << ", key=" << self.key << ", site=" << self.site
                     << ", derived_state = "
                     << static_cast<int>(self.derived_state)
                     << ", neutral = " << self.neutral << ')';
                 return out.str();
             })
        .def(py::pickle(
            [](const fwdpp::ts::mutation_record& m) {
                return py::make_tuple(m.node, m.key, m.site, m.derived_state,
                                      m.neutral);
            },
            [](py::tuple t) {
                return fwdpp::ts::mutation_record{
                    t[0].cast<decltype(fwdpp::ts::mutation_record::node)>(),
                    t[1].cast<decltype(fwdpp::ts::mutation_record::key)>(),
                    t[2].cast<decltype(fwdpp::ts::mutation_record::site)>(),
                    t[3].cast<decltype(
                        fwdpp::ts::mutation_record::derived_state)>(),
                    t[4].cast<decltype(fwdpp::ts::mutation_record::neutral)>()
                };
            }));
}

