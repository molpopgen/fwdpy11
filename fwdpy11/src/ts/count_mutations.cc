#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <fwdpy11/types/Population.hpp>
#include <fwdpp/ts/count_mutations.hpp>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<fwdpy11::Mutation>);

void
init_count_mutations(py::module& m)
{
    m.def("count_mutations",
          [](const fwdpy11::Population& pop,
             const std::vector<fwdpp::ts::TS_NODE_INT>& samples) {
              std::vector<fwdpp::uint_t>* mc
                  = new std::vector<fwdpp::uint_t>(pop.mutations.size(), 0);
              fwdpp::ts::count_mutations(pop.tables, pop.mutations, samples,
                                         *mc);
              py::capsule cap(mc, [](void* v) {
                  delete reinterpret_cast<std::vector<fwdpp::uint_t>*>(v);
              });
              return py::array(mc->size(), mc->data(), cap);
          },
          R"delim(
          Count mutation occurrences in a sample of nodes.

          :param pop: A population
          :type pop: :class:`fwdpy11.Population`
          :param samples: List of samples
          :type samples: list

          :return: Array of mutation counts
          :rtype: numpy.ndarray
          )delim");

    m.def("count_mutations",
          [](const fwdpp::ts::table_collection& tables,
             const std::vector<fwdpy11::Mutation>& mutations,
             const std::vector<fwdpp::ts::TS_NODE_INT>& samples) {
              std::vector<fwdpp::uint_t>* mc
                  = new std::vector<fwdpp::uint_t>(mutations.size(), 0);
              fwdpp::ts::count_mutations(tables, mutations, samples, *mc);
              py::capsule cap(mc, [](void* v) {
                  delete reinterpret_cast<std::vector<fwdpp::uint_t>*>(v);
              });
              return py::array(mc->size(), mc->data(), cap);
          },
          R"delim(
          Count mutation occurrences in a sample of nodes.

          :param tables: A table collection
          :type tables: :class:`fwdpy11.ts.TableCollection`
          :param mutations: Mutation list
          :type mutations: :class:`fwdpy11.VecMutation`
          :param samples: List of samples
          :type samples: list

          :return: Array of mutation counts
          :rtype: numpy.ndarray
          )delim");
}
