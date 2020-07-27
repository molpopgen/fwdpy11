#include <fwdpy11/genetic_value_to_fitness/GeneticValueIsTrait.hpp>
#include <fwdpy11/util/array_proxy.hpp>
#include <fwdpy11/numpy/array.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>

namespace py = pybind11;

class GeneticValueIsTraitTrampoline : public fwdpy11::GeneticValueIsTrait
// Trampoline class allowing custom
// GeneticValueIsTrait to be written
// in Python
{
  public:
    using fwdpy11::GeneticValueIsTrait::GeneticValueIsTrait;

    double
    operator()(const fwdpy11::DiploidGeneticValueToFitnessData input_data) const override
    {
        PYBIND11_OVERLOAD_PURE_NAME(double, fwdpy11::GeneticValueIsTrait,
                                    "__call__", operator(), input_data);
    }

    void
    update(const fwdpy11::DiploidPopulation& pop) override
    {
        PYBIND11_OVERLOAD_PURE(void, fwdpy11::GeneticValueIsTrait, update, pop);
    }

    std::shared_ptr<fwdpy11::GeneticValueToFitnessMap>
    clone() const override
    // Implementation details from pybind11 issue 1049
    {
        auto self = py::cast(this);
        auto cloned = self.attr("clone")();

        auto keep_python_state_alive = std::make_shared<py::object>(cloned);
        auto ptr = cloned.cast<GeneticValueIsTraitTrampoline*>();
        return std::shared_ptr<GeneticValueToFitnessMap>(keep_python_state_alive, ptr);
    }
};

void
init_GeneticValueIsTrait(py::module& m)
{
    py::class_<fwdpy11::GeneticValueIsTrait, fwdpy11::GeneticValueToFitnessMap,
               GeneticValueIsTraitTrampoline>(
        m, "GeneticValueIsTrait",
        "ABC for functions mapping genetic values representing traits to "
        "fitness.")
        .def(py::init<std::size_t>(), py::arg("ndim") = 1);

    py::class_<fwdpy11::DiploidGeneticValueToFitnessData>(
        m, "DiploidGeneticValueToFitnessData", py::buffer_protocol())
        .def_readonly("offspring_metadata_index",
                      &fwdpy11::DiploidGeneticValueToFitnessData::metadata_index)
        .def_property_readonly(
            "offspring_metadata",
            [](const fwdpy11::DiploidGeneticValueToFitnessData& self) {
                return py::cast<const fwdpy11::DiploidMetadata&>(
                    self.offspring_metadata.get());
            })
        .def_property_readonly(
            "parent1_metadata",
            [](const fwdpy11::DiploidGeneticValueToFitnessData& self) {
                return py::cast<const fwdpy11::DiploidMetadata&>(
                    self.parent1_metadata.get());
            })
        .def_property_readonly(
            "parent2_metadata",
            [](const fwdpy11::DiploidGeneticValueToFitnessData& self) {
                return py::cast<const fwdpy11::DiploidMetadata&>(
                    self.parent2_metadata.get());
            })
        .def_buffer([](const fwdpy11::DiploidGeneticValueToFitnessData& self) {
            return pybind11::buffer_info(
                const_cast<double*>(self.gvalues.get().data()), sizeof(double),
                pybind11::format_descriptor<double>::format(), 1,
                {self.gvalues.get().size()}, {sizeof(double)});
        });
}
