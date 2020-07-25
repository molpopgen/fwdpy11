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
#include <pybind11/stl.h>
#include <fwdpy11/genetic_values/DiploidGeneticValue.hpp>
#include <fwdpy11/genetic_value_to_fitness/GeneticValueIsTrait.hpp>
#include <fwdpy11/util/array_proxy.hpp>
#include <fwdpp/fitness_models.hpp>
#include "../genetic_value_to_fitness/GeneticValueIsTraitData.hpp"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<fwdpy11::Mutation>);
PYBIND11_MAKE_OPAQUE(std::vector<fwdpp::haploid_genome>);

struct genome_data_proxy
{
    fwdpy11::double_array_proxy effect_sizes_proxy, dominance_proxy, positions_proxy;
    fwdpy11::uint32_array_proxy smutations_proxy;
    py::object effect_sizes, dominance, smutations, positions;
    std::vector<double> effect_sizes_cpp, dominance_cpp, positions_cpp;
    py::list pymutations;
    bool filling_mutations;

    genome_data_proxy(bool fill_mutations_list)
        : effect_sizes_proxy{}, dominance_proxy{}, positions_proxy{}, smutations_proxy{},
          effect_sizes{py::cast<fwdpy11::double_array_proxy*>(&effect_sizes_proxy)},
          dominance{py::cast<fwdpy11::double_array_proxy*>(&dominance_proxy)},
          smutations{py::cast<fwdpy11::uint32_array_proxy*>(&smutations_proxy)},
          positions{py::cast<fwdpy11::double_array_proxy*>(&positions_proxy)},
          effect_sizes_cpp{}, dominance_cpp{}, positions_cpp{}, pymutations{},
          filling_mutations{fill_mutations_list}
    {
    }

    void
    set_data(const std::vector<std::uint32_t>& smutations,
             const std::vector<fwdpy11::Mutation>& mutations)
    {
        effect_sizes_cpp.clear();
        dominance_cpp.clear();
        positions_cpp.clear();
        if (filling_mutations)
            {
                pymutations.attr("clear")();
            }
        for (auto k : smutations)
            {
                auto& m = mutations[k];
                effect_sizes_cpp.push_back(m.s);
                dominance_cpp.push_back(m.h);
                positions_cpp.push_back(m.pos);
                if (filling_mutations)
                    {
                        pymutations.append(
                            py::cast<const fwdpy11::Mutation*>(&mutations[k]));
                    }
            }
        effect_sizes_proxy.set(effect_sizes_cpp);
        dominance_proxy.set(dominance_cpp);
        positions_proxy.set(positions_cpp);
        smutations_proxy.set(smutations);
    }
};

struct PyDiploidGeneticValueData
{
    py::list
    fill_list(py::object o1, py::object o2)
    {
        py::list rv;
        rv.append(o1);
        rv.append(o2);
        return rv;
    }

    fwdpy11::DiploidMetadata metadata_proxy, parent1_metadata_proxy,
        parent2_metadata_proxy;
    genome_data_proxy genome1_data, genome2_data;
    fwdpy11::double_array_proxy gvalues_proxy;
    py::object offspring_metadata, parent1_metadata, parent2_metadata, genome1, genome2,
        gvalues;
    py::list genomes, parental_metadata;
    std::size_t offspring_metadata_index;
    PyDiploidGeneticValueData(bool fill_mutations_list)
        : metadata_proxy{}, genome1_data{fill_mutations_list},
          genome2_data{fill_mutations_list},
          offspring_metadata{py::cast<fwdpy11::DiploidMetadata*>(&metadata_proxy)},
          parent1_metadata{py::cast<fwdpy11::DiploidMetadata*>(&parent1_metadata_proxy)},
          parent2_metadata{py::cast<fwdpy11::DiploidMetadata*>(&parent2_metadata_proxy)},
          genome1{py::cast<genome_data_proxy*>(&genome1_data)},
          genome2{py::cast<genome_data_proxy*>(&genome2_data)},
          gvalues{py::cast<fwdpy11::double_array_proxy*>(&gvalues_proxy)},
          genomes{fill_list(genome1, genome2)}, parental_metadata{fill_list(
                                                    parent1_metadata, parent2_metadata)},
          offspring_metadata_index{std::numeric_limits<std::size_t>::max()}
    {
    }
};

class PyDiploidGeneticValue : public fwdpy11::DiploidGeneticValue
{
  private:
    std::shared_ptr<fwdpy11::GeneticValueToFitnessMap>
    dispatch_gv2w(std::size_t ndim, py::object gvalue_to_fitness_map)
    {
        fwdpy11::GeneticValueIsFitness w(ndim);
        if (gvalue_to_fitness_map.is_none())
            {
                return w.clone();
            }
        return gvalue_to_fitness_map.cast<fwdpy11::GeneticValueIsTrait&>().clone();
    }

    std::shared_ptr<fwdpy11::GeneticValueNoise>
    dispatch_noise(py::object noise)
    {
        fwdpy11::NoNoise nonoise;
        if (noise.is_none())
            {
                return nonoise.clone();
            }
        return noise.cast<fwdpy11::GeneticValueNoise&>().clone();
    }

  public:
    bool filling_mutations;
    PyDiploidGeneticValue(std::size_t ndim, py::object gvalue_to_fitness_map,
                          py::object noise, bool fill_mutations_list)
        : fwdpy11::DiploidGeneticValue(ndim, *dispatch_gv2w(ndim, gvalue_to_fitness_map),
                                       *dispatch_noise(noise)),
          filling_mutations{fill_mutations_list}
    {
    }
};

class PyDiploidGeneticValueTrampoline : public PyDiploidGeneticValue
{
  private:
    mutable PyDiploidGeneticValueData data;
    mutable GeneticValueIsTraitData gv2w_data;
    py::object pydata;
    py::object pygv2wdata;

  public:
    PyDiploidGeneticValueTrampoline(std::size_t ndim, py::object gvalue_to_fitness_map,
                                    py::object noise, bool fill_mutations_list)
        : PyDiploidGeneticValue(ndim, gvalue_to_fitness_map, noise, fill_mutations_list),
          data{fill_mutations_list}, gv2w_data{},
          pydata{py::cast<PyDiploidGeneticValueData*>(&data)},
          pygv2wdata{py::cast<GeneticValueIsTraitData*>(&gv2w_data)}
    {
    }

    double
    calculate_gvalue(const fwdpy11::DiploidGeneticValueData input_data) override
    {
        data.metadata_proxy = input_data.offspring_metadata.get();
        data.parent1_metadata_proxy
            = input_data.pop.get()
                  .diploid_metadata[input_data.offspring_metadata.get().parents[0]];
        data.parent2_metadata_proxy
            = input_data.pop.get()
                  .diploid_metadata[input_data.offspring_metadata.get().parents[1]];
        data.genome1_data.set_data(
            input_data.pop.get()
                .haploid_genomes[input_data.pop.get()
                                     .diploids[input_data.offspring_metadata.get().label]
                                     .first]
                .smutations,
            input_data.pop.get().mutations);
        data.genome2_data.set_data(
            input_data.pop.get()
                .haploid_genomes[input_data.pop.get()
                                     .diploids[input_data.offspring_metadata.get().label]
                                     .second]
                .smutations,
            input_data.pop.get().mutations);
        data.gvalues_proxy.data = gvalues.data();
        data.gvalues_proxy.size = gvalues.size();
        data.offspring_metadata_index = input_data.metadata_index;
        PYBIND11_OVERLOAD_PURE(double, PyDiploidGeneticValue, calculate_gvalue, data);
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
                set_data(input_data, gv2w_data);
                auto obj = overload(pygv2wdata);
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
strict_additive_effects(const fwdpy11::Diploid& diploid,
                        const std::vector<fwdpp::haploid_genome>& genomes,
                        const std::vector<fwdpy11::Mutation>& mutations)
{
    double g = 0.0;
    for (auto k : genomes[diploid.first].smutations)
        {
            g += mutations[k].s;
        }
    for (auto k : genomes[diploid.second].smutations)
        {
            g += mutations[k].s;
        }
    return g;
}

double
additive_effects(const fwdpy11::Diploid& diploid,
                 const std::vector<fwdpp::haploid_genome>& genomes,
                 const std::vector<fwdpy11::Mutation>& mutations, double scaling)
{
    return fwdpp::additive_diploid(fwdpp::trait(scaling))(diploid, genomes, mutations);
}

void
init_PyDiploidGeneticValue(py::module& m)
{
    py::class_<PyDiploidGeneticValue, fwdpy11::DiploidGeneticValue,
               PyDiploidGeneticValueTrampoline>(m, "PyDiploidGeneticValue")
        .def(py::init<std::size_t, py::object, py::object, bool>(), py::arg("ndim"),
             py::arg("genetic_value_to_fitness"), py::arg("noise"),
             py::arg("fill_mutations"));

    py::class_<PyDiploidGeneticValueData>(m, "PyDiploidGeneticValueData")
        .def_readwrite("offspring_metadata",
                       &PyDiploidGeneticValueData::offspring_metadata)
        .def_readwrite("parental_metadata",
                       &PyDiploidGeneticValueData::parental_metadata)
        .def_readwrite("gvalues", &PyDiploidGeneticValueData::gvalues)
        .def_readonly("genomes", &PyDiploidGeneticValueData::genomes)
        .def_readonly("offspring_metadata_index",
                      &PyDiploidGeneticValueData::offspring_metadata_index);

    py::class_<genome_data_proxy>(m, "HaploidGenomeProxy")
        .def_readonly("effect_sizes", &genome_data_proxy::effect_sizes)
        .def_readonly("dominance", &genome_data_proxy::dominance)
        .def_readonly("positions", &genome_data_proxy::positions)
        .def_readonly("smutations", &genome_data_proxy::smutations)
        .def_property_readonly("mutations",
                               [](const genome_data_proxy& self) -> py::object {
                                   if (self.filling_mutations == false)
                                       {
                                           return py::none();
                                       }
                                   return self.pymutations;
                               });

    m.def("strict_additive_effects", &strict_additive_effects);
    m.def("additive_effects", &additive_effects);
}
