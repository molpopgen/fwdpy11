#include <fwdpp/forward_types.hpp>
#include <pybind11/stl.h>
#include <fwdpy11/numpy/array.hpp>

namespace py = pybind11;

static const auto CLASS_DOCSTRING =
    R"delim(
A haploid genome.  This object represents the ordered
mutations inherited from a single parent.
)delim";

static const auto INIT_DOCSTRING_1 =
    R"delim(
Construct haploid genome from tuple.

The tuple must be (n, mutations, smutations).
)delim";

void
init_HaploidGenome(py::module &m)
{
    py::class_<fwdpp::haploid_genome>(m, "HaploidGenome", CLASS_DOCSTRING)
        .def(py::init<fwdpp::haploid_genome::constructor_tuple>(), INIT_DOCSTRING_1)
        .def_readonly("n", &fwdpp::haploid_genome::n,
                      "Number of occurrences in the population. This has "
                      "little meaning beyond book-keeping used by the C++ "
                      "back-end. (read-only)")
        .def_property_readonly(
            "mutations",
            [](const fwdpp::haploid_genome &self) {
                return fwdpy11::make_1d_ndarray_readonly(self.mutations);
            },
            "List of keys to neutral mutations. Contains unsigned "
            "32-bit integers corresponding to mutations in the "
            "population. (read-only)")
        .def_property_readonly(
            "smutations",
            [](const fwdpp::haploid_genome &self) {
                return fwdpy11::make_1d_ndarray_readonly(self.smutations);
            },
            "List of keys to selected mutations. Contains unsigned "
            "32-bit integers corresponding to mutations in the "
            "population. (read-only)")
        .def(py::pickle(
            [](const fwdpp::haploid_genome &g) {
                return py::make_tuple(g.n, g.mutations, g.smutations);
            },
            [](py::tuple t) {
                return fwdpp::haploid_genome(t[0].cast<fwdpp::uint_t>(),
                                             t[1].cast<std::vector<fwdpp::uint_t>>(),
                                             t[2].cast<std::vector<fwdpp::uint_t>>());
            }))
        .def("__eq__", [](const fwdpp::haploid_genome &a,
                          const fwdpp::haploid_genome &b) { return a == b; });
}

