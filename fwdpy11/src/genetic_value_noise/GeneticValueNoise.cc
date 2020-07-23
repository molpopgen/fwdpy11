#include <pybind11/pybind11.h>
#include <fwdpy11/genetic_value_noise/GeneticValueNoise.hpp>

namespace py = pybind11;

struct NoiseFunctionData
{
    fwdpy11::DiploidMetadata offspring_copy, parent1_copy, parent2_copy;
    py::object offspring, parent1, parent2;
    std::size_t offspring_metadata_index;
    NoiseFunctionData()
        : offspring_copy{}, parent1_copy{}, parent2_copy{},
          offspring{py::cast<fwdpy11::DiploidMetadata*>(&offspring_copy)},
          parent1{py::cast<fwdpy11::DiploidMetadata*>(&parent1_copy)},
          parent2{py::cast<fwdpy11::DiploidMetadata*>(&parent2_copy)},
          offspring_metadata_index{std::numeric_limits<std::size_t>::max()}
    {
    }
};

class GeneticValueNoiseTrampoline : public fwdpy11::GeneticValueNoise
{
  private:
    mutable NoiseFunctionData data;
    py::object pydata;

  public:
    GeneticValueNoiseTrampoline()
        : fwdpy11::GeneticValueNoise(), data{}, pydata{
                                                    py::cast<NoiseFunctionData*>(&data)}
    {
    }

    double
    operator()(const fwdpy11::DiploidGeneticValueNoiseData input_data) const override
    {
        data.offspring_copy = input_data.offspring_metadata.get();
        data.parent1_copy = input_data.parent1_metadata.get();
        data.parent2_copy = input_data.parent2_metadata.get();
        data.offspring_metadata_index = input_data.metadata_index;
        PYBIND11_OVERLOAD_PURE_NAME(double, fwdpy11::GeneticValueNoise,
                                    "__call__", operator(), input_data.rng.get(),
                                    pydata);
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

    py::class_<NoiseFunctionData>(m, "PyGeneticValueNoiseData")
        .def_readonly("offspring_metadata", &NoiseFunctionData::offspring)
        .def_readonly("offspring_metadata_index",
                      &NoiseFunctionData::offspring_metadata_index)
        .def_readonly("parent1_metadata", &NoiseFunctionData::parent1)
        .def_readonly("parent2_metadata", &NoiseFunctionData::parent2);
}
