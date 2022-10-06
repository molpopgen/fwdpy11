#include <pybind11/pybind11.h>
#include <fwdpy11/genetic_value_noise/GeneticValueNoise.hpp>

namespace py = pybind11;

class GeneticValueNoiseTrampoline : public fwdpy11::GeneticValueNoise
{
  public:
    using fwdpy11::GeneticValueNoise::GeneticValueNoise;

    double
    operator()(const fwdpy11::DiploidGeneticValueNoiseData input_data) const override
    {
        PYBIND11_OVERLOAD_PURE_NAME(double, fwdpy11::GeneticValueNoise,
                                    "__call__", operator(), input_data);
    }

    void
    update(const fwdpy11::DiploidPopulation& pop) override
    {
        PYBIND11_OVERLOAD_PURE(void, fwdpy11::GeneticValueNoise, update, pop);
    }

    std::shared_ptr<fwdpy11::GeneticValueNoise>
    clone() const override
    // Implementation details from pybind11 issue 1049
    {
        auto self = py::cast(this);
        auto cloned = self.attr("clone")();

        auto keep_python_state_alive = std::make_shared<py::object>(cloned);
        auto ptr = cloned.cast<GeneticValueNoiseTrampoline*>();
        return std::shared_ptr<GeneticValueNoise>(keep_python_state_alive, ptr);
    }
};

void
init_GeneticValueNoise(py::module& m)
{
    py::class_<fwdpy11::GeneticValueNoise, GeneticValueNoiseTrampoline>(
        m, "GeneticValueNoise",
        "ABC for noise classes affecting :class:`fwdpy11.DiploidPopulation`.")
        .def(py::init<>());

    py::class_<fwdpy11::DiploidGeneticValueNoiseData>(m, "DiploidGeneticValueNoiseData")
        .def_property_readonly("rng",
                               [](const fwdpy11::DiploidGeneticValueNoiseData& self) {
                                   return py::cast<const fwdpy11::GSLrng_t&>(
                                       self.rng.get());
                               })
        .def_readonly("offspring_metadata_index",
                      &fwdpy11::DiploidGeneticValueNoiseData::metadata_index)
        .def_property_readonly("offspring_metadata",
                               [](const fwdpy11::DiploidGeneticValueNoiseData& self) {
                                   return py::cast<const fwdpy11::DiploidMetadata&>(
                                       self.offspring_metadata.get());
                               })
        .def_property_readonly("parent1_metadata",
                               [](const fwdpy11::DiploidGeneticValueNoiseData& self) {
                                   return py::cast<const fwdpy11::DiploidMetadata&>(
                                       self.parent1_metadata.get());
                               })
        .def_property_readonly("parent2_metadata",
                               [](const fwdpy11::DiploidGeneticValueNoiseData& self) {
                                   return py::cast<const fwdpy11::DiploidMetadata&>(
                                       self.parent2_metadata.get());
                               });
}
