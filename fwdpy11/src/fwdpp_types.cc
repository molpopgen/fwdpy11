#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/sugar/popgenmut.hpp>

namespace py = pybind11;

PYBIND11_PLUGIN(fwdpp_types) {
    py::module m("fwdpp_types", "example extending");

    // low-level types

    py::class_<KTfwd::mutation_base>(m, "MutationBase",
R"delim(
Base class for mutations.
)delim")
        .def(py::init<double, bool, std::uint16_t>(), "Constructor")
        .def_readwrite("pos", &KTfwd::mutation_base::pos, "Position (float).")
        .def_readwrite("neutral", &KTfwd::mutation_base::neutral, "Boolean")
        .def_readwrite("xtra", &KTfwd::mutation_base::xtra,
                       "16-bits worth of extra data.");

    py::class_<KTfwd::gamete>(m, "Gamete")
        .def_readonly("n", &KTfwd::gamete::n,"Number of occurrences in the population.")
        .def_readonly("mutations", &KTfwd::gamete::mutations,"Vector of keys to neutral mutations.")
        .def_readonly("smutations", &KTfwd::gamete::smutations,"Vector of keys to selected mutations.")
        .def("as_dict",
             // This lambda shows that
             // we can return dicts
             // with a mix of scalars + containers
             [](const KTfwd::gamete &g) noexcept {
        using obj = pybind11::object;
        pybind11::dict rv;
        rv[obj(pybind11::cast("n"))] = obj(pybind11::cast(g.n));
        rv[obj(pybind11::cast("mutations"))] = obj(pybind11::cast(g.mutations));
        rv[obj(pybind11::cast("smutations"))] =
            obj(pybind11::cast(g.smutations));
        return rv;
             })
        .def("__getstate__",[](const KTfwd::gamete &g)
        {
        return py::make_tuple(g.n, g.mutations, g.smutations);
        })
        .def("__setstate__",[](KTfwd::gamete & g, py::tuple t)
        {
        new (&g) KTfwd::gamete(t[0].cast<KTfwd::uint_t>(),
                               t[1].cast<std::vector<KTfwd::uint_t>>(),
                               t[2].cast<std::vector<KTfwd::uint_t>>());
        });
    
    // Sugar types
    py::class_<KTfwd::popgenmut, KTfwd::mutation_base>(
        m, "Mutation", "Mutation with effect size and dominance")
        .def(py::init<double, double, double, unsigned, std::uint16_t>())
        .def_readwrite("g", &KTfwd::popgenmut::g,"Generation when mutation arose (origination time).")
        .def_readwrite("s", &KTfwd::popgenmut::s,"Selection coefficient/effect size.")
        .def_readwrite("h", &KTfwd::popgenmut::h,"Dominance/effect in heterozygotes.")
        .def("__getstate__",
             [](const KTfwd::popgenmut &m) {
                 return py::make_tuple(m.pos, m.s, m.h, m.g, m.xtra);
             })
        .def("__setstate__", [](KTfwd::popgenmut &m, py::tuple p) {
            new (&m) KTfwd::popgenmut(
                p[0].cast<double>(), p[1].cast<double>(), p[2].cast<double>(),
                p[3].cast<unsigned>(), p[4].cast<std::uint16_t>());
        });

    return m.ptr();
}
