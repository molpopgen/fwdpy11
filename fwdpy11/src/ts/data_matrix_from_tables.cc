#include <fwdpp/ts/generate_data_matrix.hpp>
#include <fwdpy11/types/Population.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void
init_data_matrix_from_tables(py::module& m)
{
    m.def("make_data_matrix",
          [](const fwdpy11::Population& pop,
             const std::vector<fwdpp::ts::TS_NODE_INT>& samples,
             const bool record_neutral, const bool record_selected) {
              return fwdpp::ts::generate_data_matrix(
                  pop.tables, samples, pop.mutations, record_neutral,
                  record_selected);
          },
          py::arg("pop"), py::arg("samples"), py::arg("record_neutral"),
          py::arg("record_selected"),
          R"delim(
     Create a :class:`fwdpy11.sampling.DataMatrix` from a table collection.
     
     :param pop: A population
     :type pop: :class:`fwdpy11.Population`
     :param samples: A list of sample nodes
     :type samples: list
     :param record_neutral: If True, generate data for neutral variants
     :type record_neutral: boolean
     :param record_selected: If True, generate data for selected variants
     :type record_selected: boolean

     .. deprecated:: 0.3.0

        Prefer :func:`fwdpy11.ts.data_matrix_from_tables`.
     )delim");

    m.def("data_matrix_from_tables",
          [](const fwdpp::ts::table_collection& tables,
             const std::vector<fwdpy11::Mutation>& mutations,
             const std::vector<fwdpp::ts::TS_NODE_INT>& samples,
             bool record_neutral, bool record_selected) {
              return fwdpp::ts::generate_data_matrix(
                  tables, samples, mutations, record_neutral, record_selected);
          },
          R"delim(
     Create a :class:`fwdpy11.sampling.DataMatrix` from a table collection.
     
     :param tables: A TableCollection
     :type tables: :class:`fwdpy11.ts.TableCollection`
     :param mutations: Container of mutations
     :type mutations: :class:`fwdpy11.VecMutation`
     :param samples: A list of sample nodes
     :type samples: list or array
     :param record_neutral: If True, generate data for neutral variants
     :type record_neutral: boolean
     :param record_selected: If True, generate data for selected variants
     :type record_selected: boolean

     :rtype: :class:`fwdpy11.sampling.DataMatrix`

     .. versionadded:: 0.3.0
     )delim");

}
