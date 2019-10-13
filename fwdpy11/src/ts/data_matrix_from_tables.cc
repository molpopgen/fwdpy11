#include <fwdpp/ts/generate_data_matrix.hpp>
#include <fwdpy11/types/Population.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<fwdpy11::Mutation>);

fwdpp::data_matrix
generate_data_matrix(const fwdpp::ts::table_collection& tables,
                     const std::vector<fwdpp::ts::TS_NODE_INT>& samples,
                     bool record_neutral, bool record_selected, bool include_fixations, double start,
                     double stop)
{
    return fwdpp::ts::generate_data_matrix(tables, samples,
                                           record_neutral, record_selected, !include_fixations,
                                           start, stop);
}

void
init_data_matrix_from_tables(py::module& m)
{
    m.def(
        "make_data_matrix",
        [](const fwdpy11::Population& pop,
           const std::vector<fwdpp::ts::TS_NODE_INT>& samples,
           const bool record_neutral, const bool record_selected) {
            return fwdpp::ts::generate_data_matrix(
                pop.tables, samples, record_neutral,
                record_selected, true);
        },
        py::arg("pop"), py::arg("samples"), py::arg("record_neutral"),
        py::arg("record_selected"),
        R"delim(
     Create a :class:`fwdpy11.DataMatrix` from a table collection.
     
     :param pop: A population
     :type pop: :class:`fwdpy11.Population`
     :param samples: A list of sample nodes
     :type samples: list
     :param record_neutral: If True, generate data for neutral variants
     :type record_neutral: boolean
     :param record_selected: If True, generate data for selected variants
     :type record_selected: boolean

     .. deprecated:: 0.3.0

        Prefer :func:`fwdpy11.data_matrix_from_tables`.
     )delim");

    m.def("data_matrix_from_tables", &generate_data_matrix, py::arg("tables"),
          py::arg("samples"), py::arg("record_neutral"),
          py::arg("record_selected"),
          py::arg("include_fixations") = false, py::arg("begin") = 0.0,
          py::arg("end") = std::numeric_limits<double>::max(),
          R"delim(
     Create a :class:`fwdpy11.DataMatrix` from a table collection.
     
     :param tables: A TableCollection
     :type tables: :class:`fwdpy11.TableCollection`
     :param samples: A list of sample nodes
     :type samples: list or array
     :param record_neutral: If True, generate data for neutral variants
     :type record_neutral: boolean
     :param record_selected: If True, generate data for selected variants
     :type record_selected: boolean
     :param include_selected: (False) Whether to include variants fixed in the sample
     :type include_selected: boolean
     :param begin: (0.0) Start of range, inclusive
     :param end: (max float) End of range, exclusive

     :rtype: :class:`fwdpy11.DataMatrix`

     .. versionadded:: 0.3.0

     .. versionchanged:: 0.4.1
        
            Add begin, end options as floats

     .. versionchanged:: 0.5.0

            No longer requires :class:`fwdpy11.MutationVector` argument
     )delim");
}
