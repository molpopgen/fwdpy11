#include <cstdint>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <core/demes/forward_graph.hpp>

namespace py = pybind11;

void
init_forward_graph(py::module &m)
{
    py::class_<fwdpy11_core::ForwardDemesGraph>(m, "_ForwardDemesGraph")
        .def(py::init<const std::string & /*yaml*/, std::uint32_t /*burnin*/,
                      bool /*round_epoch_sizes*/>(),
             py::arg("yaml"), py::arg("burnin"), py::arg("round_epoch_sizes"),
             py::kw_only())
        .def("_sum_deme_sizes_at_time_zero",
             &fwdpy11_core::ForwardDemesGraph::sum_deme_sizes_at_time_zero)
        .def("_model_end_time", &fwdpy11_core::ForwardDemesGraph::model_end_time)
        .def("_parental_deme_sizes_at_time_zero",
             &fwdpy11_core::ForwardDemesGraph::parental_deme_sizes_at_time_zero);
}
