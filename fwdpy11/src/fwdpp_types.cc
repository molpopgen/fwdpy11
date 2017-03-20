#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/sugar/popgenmut.hpp>

namespace py = pybind11;

PYBIND11_PLUGIN(fwdpp_types) {
    py::module m("fwdpp_types", "example extending");

    // low-level types

    py::class_<KTfwd::mutation_base>(m, "MutationBase",
                                     "Base class for mutations.")
        .def(py::init<double, bool, std::uint16_t>(), "Constructor")
        .def_readwrite("pos", &KTfwd::mutation_base::pos, "Position (float).")
        .def_readwrite("neutral", &KTfwd::mutation_base::neutral, "Boolean")
        .def_readwrite("xtra", &KTfwd::mutation_base::xtra,
                       "16-bits worth of extra data");

    py::class_<KTfwd::gamete>(m, "Gamete")
        .def_readonly("n", &KTfwd::gamete::n)
        .def_readonly("mutations", &KTfwd::gamete::mutations)
        .def_readonly("smutations", &KTfwd::gamete::smutations)
        .def("as_dict",
             // This lambda shows that
             // we can return dicts
             // with a mix of scalars + containers
             [](const KTfwd::gamete &g) noexcept {
                 using obj = pybind11::object;
                 pybind11::dict rv;
                 rv[obj(pybind11::cast("n"))] = obj(pybind11::cast(g.n));
                 rv[obj(pybind11::cast("mutations"))] =
                     obj(pybind11::cast(g.mutations));
                 rv[obj(pybind11::cast("smutations"))] =
                     obj(pybind11::cast(g.smutations));
                 return rv;
             });
    ;

    // Sugar types
    py::class_<KTfwd::popgenmut, KTfwd::mutation_base>(
        m, "Mutation", "Mutation with effect size and dominance")
        .def(py::init<double, double, double, unsigned, std::uint16_t>())
        .def_readwrite("g", &KTfwd::popgenmut::g)
        .def_readwrite("s", &KTfwd::popgenmut::s)
        .def_readwrite("h", &KTfwd::popgenmut::h);

    return m.ptr();
}
