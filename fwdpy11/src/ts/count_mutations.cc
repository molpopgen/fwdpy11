#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <fwdpy11/types/Population.hpp>
#include <fwdpy11/numpy/array.hpp>
#include <fwdpp/ts/count_mutations.hpp>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<fwdpy11::Mutation>);

void
init_count_mutations(py::module& m)
{
    m.def(
        "count_mutations",
        [](const fwdpy11::Population& pop,
           const std::vector<fwdpp::ts::table_index_t>& samples) {
            std::vector<fwdpp::uint_t> mc(pop.mutations.size(), 0);
            fwdpp::ts::count_mutations(*pop.tables, pop.mutations, samples, mc);
            return fwdpy11::make_1d_array_with_capsule(std::move(mc));
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

    m.def(
        "count_mutations",
        [](const fwdpp::ts::std_table_collection& tables,
           const std::vector<fwdpy11::Mutation>& mutations,
           const std::vector<fwdpp::ts::table_index_t>& samples) {
            std::vector<fwdpp::uint_t> mc(mutations.size(), 0);
            fwdpp::ts::count_mutations(tables, mutations, samples, mc);
            return fwdpy11::make_1d_array_with_capsule(std::move(mc));
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
