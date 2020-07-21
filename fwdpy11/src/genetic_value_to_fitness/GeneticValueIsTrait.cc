#include <fwdpy11/genetic_value_to_fitness/GeneticValueIsTrait.hpp>
#include <fwdpy11/util/array_proxy.hpp>
#include <fwdpy11/numpy/array.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include "GeneticValueIsTraitData.hpp"

namespace py = pybind11;


class GeneticValueIsTraitTrampoline : public fwdpy11::GeneticValueIsTrait
// Trampoline class allowing custom
// GeneticValueIsTrait to be written
// in Python
{
  private:
    mutable GeneticValueIsTraitData data;
    py::object pydata;

  public:
    GeneticValueIsTraitTrampoline(std::size_t ndim)
        : fwdpy11::GeneticValueIsTrait(ndim), data{},
          pydata{py::cast<GeneticValueIsTraitData*>(&data)}
    {
    }

    double
    operator()(const fwdpy11::DiploidGeneticValueToFitnessData input_data) const override
    {
        data.offspring_metadata_copy = input_data.offspring_metadata.get();
        data.buffer.data = const_cast<double*>(input_data.gvalues.get().data());
        data.buffer.size = input_data.gvalues.get().size();
        PYBIND11_OVERLOAD_PURE_NAME(double, fwdpy11::GeneticValueIsTrait,
                                    "__call__", operator(), pydata);
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

    py::class_<GeneticValueIsTraitData>(m, "PyGeneticValueIsTraitData")
        .def_readonly("offspring_metadata", &GeneticValueIsTraitData::offspring_metadata)
        .def_readonly("genetic_values", &GeneticValueIsTraitData::genetic_values);
}
