#include <cstdint>
#include <pybind11/pybind11.h>
#include <core/demes/forward_graph.hpp>

namespace py = pybind11;

void
init_forward_graph(py::module &m)
{
    py::class_<fwdpy11_core::ForwardDemesGraph>(m, "_ForwardDemesGraph")
        .def(py::init<const std::string & /*yaml*/, std::uint32_t /*burnin*/,
                      bool /*round_epoch_sizes*/>(),
             py::arg("yaml"), py::arg("burnin"), py::arg("round_epoch_sizes"),
             py::kw_only());
}
