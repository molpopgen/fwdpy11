#include <pybind11/pybind11.h>
#include <fwdpy11/types/Diploid.hpp>
namespace py = pybind11;

void
init_DiploidGenotype(py::module &m)
{
    py::class_<fwdpy11::DiploidGenotype>(
        m, "DiploidGenotype",
        "Diploid data type for a single (usually contiguous) genomic region")
        .def(py::init<>())
        .def(py::init<std::size_t, std::size_t>())
        .def_readonly("first", &fwdpy11::DiploidGenotype::first,
                      "Key to first gamete. (read-only)")
        .def_readonly("second", &fwdpy11::DiploidGenotype::second,
                      "Key to second gamete. (read-only)")
        .def(py::pickle(
            [](const fwdpy11::DiploidGenotype &d) {
                return py::make_tuple(d.first, d.second);
            },
            [](py::tuple t) {
                std::unique_ptr<fwdpy11::DiploidGenotype> d(
                    new fwdpy11::DiploidGenotype{ t[0].cast<std::size_t>(),
                                                  t[1].cast<std::size_t>() });
                return d;
            }))
        .def("__eq__",
             [](const fwdpy11::DiploidGenotype &a,
                const fwdpy11::DiploidGenotype &b) { return a == b; });
}
