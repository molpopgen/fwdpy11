#include <pybind11/pybind11.h>
#include <fwdpy11/genetic_value_noise/GeneticValueNoise.hpp>

namespace py = pybind11;

class GeneticValueNoiseTrampoline : public fwdpy11::GeneticValueNoise
{
  public:
    using fwdpy11::GeneticValueNoise::GeneticValueNoise;

    double
    operator()(const fwdpy11::GSLrng_t& rng,
               const fwdpy11::DiploidMetadata& offspring_metadata,
               const std::size_t parent1, const std::size_t parent2,
               const fwdpy11::DiploidPopulation& pop) const override
    {
        PYBIND11_OVERLOAD_PURE_NAME(double, fwdpy11::GeneticValueNoise,
                                    "__call__", operator(), rng, offspring_metadata, parent1,
                                    parent2, pop);
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
}
