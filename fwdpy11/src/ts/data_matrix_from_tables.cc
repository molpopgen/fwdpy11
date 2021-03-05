#include <fwdpp/ts/generate_data_matrix.hpp>
#include <fwdpy11/types/Population.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<fwdpy11::Mutation>);

fwdpp::data_matrix
generate_data_matrix(const fwdpp::ts::std_table_collection& tables,
                     const std::vector<fwdpp::ts::table_index_t>& samples,
                     bool record_neutral, bool record_selected, bool include_fixations,
                     double start, double stop)
{
    return fwdpp::ts::generate_data_matrix(tables, samples, record_neutral,
                                           record_selected, !include_fixations, start,
                                           stop);
}

void
init_data_matrix_from_tables(py::module& m)
{
    m.def(
        "_make_data_matrix",
        [](const fwdpy11::Population& pop,
           const std::vector<fwdpp::ts::table_index_t>& samples,
           const bool record_neutral, const bool record_selected) {
            return fwdpp::ts::generate_data_matrix(*pop.tables, samples, record_neutral,
                                                   record_selected, true);
        },
        py::arg("pop"), py::arg("samples"), py::arg("record_neutral"),
        py::arg("record_selected"));

    m.def("_data_matrix_from_tables", &generate_data_matrix, py::arg("tables"),
          py::arg("samples"), py::arg("record_neutral"), py::arg("record_selected"),
          py::arg("include_fixations") = false, py::arg("begin") = 0.0,
          py::arg("end") = std::numeric_limits<double>::max());
}
