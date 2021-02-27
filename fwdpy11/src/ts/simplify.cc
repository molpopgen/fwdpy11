#include <fwdpy11/types/Population.hpp>
#include <fwdpy11/numpy/array.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <fwdpp/ts/table_simplifier.hpp>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<fwdpy11::Mutation>);

py::tuple
simplify(const fwdpy11::Population& pop,
         const std::vector<fwdpp::ts::table_index_t>& samples)
{
    if (pop.tables->genome_length() == std::numeric_limits<double>::max())
        {
            throw std::invalid_argument("population is not using tree sequences");
        }
    if (pop.tables->num_nodes() == 0)
        {
            throw std::invalid_argument("population has empty TableCollection");
        }
    if (samples.empty())
        {
            throw std::invalid_argument("empty sample list");
        }
    if (std::any_of(samples.begin(), samples.end(),
                    [&pop](const fwdpp::ts::table_index_t s) {
                        return s == fwdpp::ts::NULL_INDEX
                               || static_cast<std::size_t>(s) >= pop.tables->num_nodes();
                    }))
        {
            throw std::invalid_argument("invalid sample list");
        }
    auto t(*pop.tables);
    fwdpp::ts::table_simplifier<fwdpp::ts::std_table_collection> simplifier{};
    auto rv = simplifier.simplify(t, samples);
    t.build_indexes();
    return py::make_tuple(std::move(t),
                          fwdpy11::make_1d_array_with_capsule(std::move(rv.first)));
}

void
init_simplify_functions(py::module& m)
{
    m.def("_simplify", &simplify, py::arg("pop"), py::arg("samples"));

    m.def(
        "_simplify_tables",
        [](const fwdpp::ts::std_table_collection& tables,
           const std::vector<fwdpp::ts::table_index_t>& samples) -> py::tuple {
            auto t(tables);
            fwdpp::ts::table_simplifier<fwdpp::ts::std_table_collection> simplifier{};
            auto rv = simplifier.simplify(t, samples);
            t.build_indexes();
            return py::make_tuple(
                std::move(t), fwdpy11::make_1d_array_with_capsule(std::move(rv.first)));
        },
        py::arg("tables"), py::arg("samples"));
}

