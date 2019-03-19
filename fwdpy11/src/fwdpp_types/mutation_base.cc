#include <fwdpp/forward_types.hpp>
#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_mutation_base(py::module & m)
{
    py::class_<fwdpp::mutation_base>(m, "MutationBase",
                                     R"delim(
                                        Base class for mutations.
                                     )delim")
        .def(py::init<double, bool, std::uint16_t>(), "Constructor")
        .def_readonly("pos", &fwdpp::mutation_base::pos, "Position (float).")
        .def_readonly("neutral", &fwdpp::mutation_base::neutral, "Boolean")
        .def_readwrite("label", &fwdpp::mutation_base::xtra,
                       "A 16-bit unsigned integer that can be used for adding "
                       "\"meta-data\" to mutations");
}
