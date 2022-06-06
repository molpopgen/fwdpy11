#include <pybind11/pybind11.h>
#include <fwdpy11/discrete_demography/forward_demes_graph.hpp>

namespace py = pybind11;

using namespace fwdpy11::discrete_demography;

void
init_ForwardDemesGraph(py::module& m)
{
    py::class_<ForwardDemesGraph>(m, "_ll_ForwardDemesGraph")
        .def(py::init<>())
        .def("_add_deme", &::ForwardDemesGraph::add_deme);

    py::class_<Selfing>(m, "_ll_Selfing");
    py::class_<SizeFunction>(m, "_ll_SizeFunction");

    py::class_<DemeRef>(m, "_ll_DemeRef")
        .def(
            "_add_epoch",
            [](DemeRef& self, demes_model_time end_time, std::uint32_t start_size,
               std::uint32_t end_size, double cloning_rate, Selfing selfing,
               SizeFunction size_function) {
                return self.deme.get().add_epoch(end_time, start_size, end_size,
                                                 cloning_rate, selfing, size_function);
            },
            py::kw_only(), py::arg("end_time"), py::arg("start_size"),
            py::arg("end_size"), py::arg("cloning_rate"), py::arg("selfing"),
            py::arg("size_function"));

    m.def("_constant_size_function", &constant_size_function);

    m.def("_wright_fisher_selfing", &Selfing::wright_fisher);
}
