//
// Copyright (C) 2017-2020 Kevin Thornton <krthornt@uci.edu>
//
// This file is part of fwdpy11.
//
// fwdpy11 is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// fwdpy11 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with fwdpy11.  If not, see <http://www.gnu.org/licenses/>.
//
#include <pybind11/pybind11.h>
#include <fwdpy11/genetic_values/DiploidGeneticValue.hpp>
#include <fwdpy11/genetic_value_to_fitness/GeneticValueIsTrait.hpp>
#include <fwdpp/fitness_models.hpp>

namespace py = pybind11;

struct PyDiploidGeneticValueData
{
    std::reference_wrapper<const fwdpy11::GSLrng_t> rng;
    std::reference_wrapper<const fwdpy11::DiploidPopulation> pop;
    std::reference_wrapper<const fwdpy11::DiploidMetadata> offspring_metadata,
        parent1_metadata, parent2_metadata;
    std::size_t metadata_index;
    std::reference_wrapper<std::vector<double>> gvalues;

    PyDiploidGeneticValueData(const fwdpy11::DiploidGeneticValueData& input_data,
                              std::vector<double>& gv)
        : rng{input_data.rng}, pop{input_data.pop},
          offspring_metadata{input_data.offspring_metadata},
          parent1_metadata{input_data.parent1_metadata},
          parent2_metadata{input_data.parent2_metadata},
          metadata_index{input_data.metadata_index}, gvalues{gv}
    {
    }
};

class PyDiploidGeneticValue : public fwdpy11::DiploidGeneticValue
{
  public:
    PyDiploidGeneticValue(std::size_t ndim,
                          const fwdpy11::GeneticValueToFitnessMap* gvalue_to_fitness_map,
                          const fwdpy11::GeneticValueNoise* noise)
        : fwdpy11::DiploidGeneticValue(ndim, gvalue_to_fitness_map, noise)
    {
    }
};

class PyDiploidGeneticValueTrampoline : public PyDiploidGeneticValue
{
  public:
    using PyDiploidGeneticValue::PyDiploidGeneticValue;

    double
    calculate_gvalue(const fwdpy11::DiploidGeneticValueData input_data) override
    {
        PYBIND11_OVERLOAD_PURE(double, PyDiploidGeneticValue, calculate_gvalue,
                               PyDiploidGeneticValueData(input_data, this->gvalues));
    }

    double
    genetic_value_to_fitness(
        const fwdpy11::DiploidGeneticValueToFitnessData input_data) override
    // NOTE: see https://pybind11.readthedocs.io/en/stable/advanced/classes.html#extended-trampoline-class-functionality
    {
        pybind11::gil_scoped_acquire gil; // Acquire the GIL while in this scope.
        // Try to look up the overloaded method on the Python side.
        pybind11::function overload
            = pybind11::get_overload(this, "genetic_value_to_fitness");
        if (overload)
            {
                auto obj = overload(input_data);
                return obj.cast<double>();
            }
        return this->gv2w->operator()(input_data);
    }

    void
    update(const fwdpy11::DiploidPopulation& pop) override
    {
        PYBIND11_OVERLOAD_PURE(void, PyDiploidGeneticValue, update, pop);
    }
};

double
strict_additive_effects(const fwdpy11::DiploidPopulation& pop,
                        const fwdpy11::DiploidMetadata& individual)
{
    double g = 0.0;
    const auto& diploid = pop.diploids[individual.label];
    for (auto k : pop.haploid_genomes[diploid.first].smutations)
        {
            g += pop.mutations[k].s;
        }
    for (auto k : pop.haploid_genomes[diploid.second].smutations)
        {
            g += pop.mutations[k].s;
        }
    return g;
}

double
additive_effects(const fwdpy11::DiploidPopulation& pop,
                 const fwdpy11::DiploidMetadata& individual, double scaling)
{
    return fwdpp::additive_diploid(fwdpp::trait(scaling))(
        pop.diploids[individual.label], pop.haploid_genomes, pop.mutations);
}

void
init_PyDiploidGeneticValue(py::module& m)
{
    py::class_<PyDiploidGeneticValue, fwdpy11::DiploidGeneticValue,
               PyDiploidGeneticValueTrampoline>(m, "PyDiploidGeneticValue")
        .def(py::init<std::size_t, const fwdpy11::GeneticValueToFitnessMap*,
                      const fwdpy11::GeneticValueNoise*>(),
             py::arg("ndim"), py::arg("genetic_value_to_fitness").none(true),
             py::arg("noise").none(true));

    py::class_<PyDiploidGeneticValueData>(m, "PyDiploidGeneticValueData",
                                          py::buffer_protocol())
        .def_property_readonly("rng",
                               [](const PyDiploidGeneticValueData& self) {
                                   return py::cast<const fwdpy11::GSLrng_t&>(
                                       self.rng.get());
                               })
        .def_property_readonly("pop",
                               [](const PyDiploidGeneticValueData& self) {
                                   return py::cast<const fwdpy11::DiploidPopulation&>(
                                       self.pop.get());
                               })
        .def_readonly("offspring_metadata_index",
                      &PyDiploidGeneticValueData::metadata_index)
        .def_property_readonly("offspring_metadata",
                               [](const PyDiploidGeneticValueData& self) {
                                   return py::cast<const fwdpy11::DiploidMetadata&>(
                                       self.offspring_metadata.get());
                               })
        .def_property_readonly("parent1_metadata",
                               [](const PyDiploidGeneticValueData& self) {
                                   return py::cast<const fwdpy11::DiploidMetadata&>(
                                       self.parent1_metadata.get());
                               })
        .def_property_readonly("parent2_metadata",
                               [](const PyDiploidGeneticValueData& self) {
                                   return py::cast<const fwdpy11::DiploidMetadata&>(
                                       self.parent2_metadata.get());
                               })
        .def_buffer([](const PyDiploidGeneticValueData& self) {
            return pybind11::buffer_info(
                const_cast<double*>(self.gvalues.get().data()), sizeof(double),
                pybind11::format_descriptor<double>::format(), 1,
                {self.gvalues.get().size()}, {sizeof(double)});
        });

    m.def("strict_additive_effects", &strict_additive_effects);
    m.def("additive_effects", &additive_effects);
}
