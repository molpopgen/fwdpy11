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

namespace py = pybind11;

template <typename T> struct array_proxy
{
    T* data;
    std::size_t size;
    array_proxy() : data{nullptr}, size{0}
    {
    }

    void
    set(const std::vector<T>& v)
    {
        data = const_cast<T*>(v.data());
        size = v.size();
    }
};

using double_array_proxy = array_proxy<double>;
using mutation_key_array_proxy = array_proxy<std::uint32_t>;

struct genome_data_proxy
{
    double_array_proxy effect_sizes_proxy, dominance_proxy, positions_proxy;
    mutation_key_array_proxy smutations_proxy;
    py::object effect_sizes, dominance, smutations, positions;
    std::vector<double> effect_sizes_cpp, dominance_cpp, positions_cpp;

    genome_data_proxy()
        : effect_sizes_proxy{}, dominance_proxy{}, positions_proxy{}, smutations_proxy{},
          effect_sizes{py::cast<double_array_proxy*>(&effect_sizes_proxy)},
          dominance{py::cast<double_array_proxy*>(&dominance_proxy)},
          smutations{py::cast<mutation_key_array_proxy*>(&smutations_proxy)},
          positions{py::cast<double_array_proxy*>(&positions_proxy)}, effect_sizes_cpp{},
          dominance_cpp{}, positions_cpp{}
    {
    }

    void
    set_data(const std::vector<std::uint32_t>& smutations,
             const std::vector<fwdpy11::Mutation>& mutations)
    {
        effect_sizes_cpp.clear();
        dominance_cpp.clear();
        positions_cpp.clear();
        for (auto k : smutations)
            {
                auto& m = mutations[k];
                effect_sizes_cpp.push_back(m.s);
                dominance_cpp.push_back(m.h);
                positions_cpp.push_back(m.pos);
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
    double_array_proxy gvalues_proxy;
    py::object offspring_metadata, parent1_metadata, parent2_metadata, genome1, genome2,
        gvalues;
    py::list genomes, parental_metadata;
    PyDiploidGeneticValueData()
        : metadata_proxy{}, genome1_data{}, genome2_data{},
          offspring_metadata{py::cast<fwdpy11::DiploidMetadata*>(&metadata_proxy)},
          parent1_metadata{py::cast<fwdpy11::DiploidMetadata*>(&parent1_metadata_proxy)},
          parent2_metadata{py::cast<fwdpy11::DiploidMetadata*>(&parent2_metadata_proxy)},
          genome1{py::cast<genome_data_proxy*>(&genome1_data)},
          genome2{py::cast<genome_data_proxy*>(&genome2_data)},
          gvalues{py::cast<double_array_proxy*>(&gvalues_proxy)}, genomes{fill_list(
                                                                      genome1, genome2)},
          parental_metadata{fill_list(parent1_metadata, parent2_metadata)}
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
    PyDiploidGeneticValue(std::size_t ndim, py::object gvalue_to_fitness_map,
                          py::object noise)
        : fwdpy11::DiploidGeneticValue(ndim, *dispatch_gv2w(ndim, gvalue_to_fitness_map),
                                       *dispatch_noise(noise))
    {
    }
};

class PyDiploidGeneticValueTrampoline : public PyDiploidGeneticValue
{
  private:
    mutable PyDiploidGeneticValueData data;
    py::object pydata;

  public:
    PyDiploidGeneticValueTrampoline(std::size_t ndim, py::object gvalue_to_fitness_map,
                                    py::object noise)
        : PyDiploidGeneticValue(ndim, gvalue_to_fitness_map, noise), data{},
          pydata{py::cast<PyDiploidGeneticValueData*>(&data)}
    {
    }

    double
    calculate_gvalue(const std::size_t diploid_index,
                     const fwdpy11::DiploidMetadata& diploid_metadata,
                     const fwdpy11::DiploidPopulation& pop) const override
    {
        data.metadata_proxy = diploid_metadata;
        data.parent1_metadata_proxy = pop.diploid_metadata[diploid_metadata.parents[0]];
        data.parent2_metadata_proxy = pop.diploid_metadata[diploid_metadata.parents[1]];
        data.genome1_data.set_data(
            pop.haploid_genomes[pop.diploids[diploid_index].first].smutations,
            pop.mutations);
        data.genome2_data.set_data(
            pop.haploid_genomes[pop.diploids[diploid_index].second].smutations,
            pop.mutations);
        data.gvalues_proxy.data = gvalues.data();
        data.gvalues_proxy.size = gvalues.size();
        PYBIND11_OVERLOAD_PURE(double, PyDiploidGeneticValue, calculate_gvalue, data);
    }

    double
    genetic_value_to_fitness(
        const fwdpy11::DiploidMetadata& offspring_metadata) const override
    {
        PYBIND11_OVERLOAD(double, PyDiploidGeneticValue, genetic_value_to_fitness,
                          offspring_metadata);
    }

    void
    update(const fwdpy11::DiploidPopulation& pop) override
    {
        PYBIND11_OVERLOAD(void, PyDiploidGeneticValue, update, pop);
        gv2w->update(pop);
        noise_fxn->update(pop);
    }
};

void
init_PyDiploidGeneticValue(py::module& m)
{
    py::class_<PyDiploidGeneticValue, fwdpy11::DiploidGeneticValue,
               PyDiploidGeneticValueTrampoline>(m, "PyDiploidGeneticValue")
        .def(py::init<std::size_t, py::object, py::object>())
        .def("update_members",
             [](PyDiploidGeneticValue& self, const fwdpy11::DiploidPopulation& pop) {
                 self.gv2w->update(pop);
                 self.noise_fxn->update(pop);
             });

    py::class_<PyDiploidGeneticValueData>(m, "PyDiploidGeneticValueData")
        .def_readwrite("offspring_metadata",
                       &PyDiploidGeneticValueData::offspring_metadata)
        .def_readwrite("parent1_metadata",
                       &PyDiploidGeneticValueData::parent1_metadata)
        .def_readwrite("parental_metadata",
                       &PyDiploidGeneticValueData::parental_metadata)
        .def_readwrite("gvalues", &PyDiploidGeneticValueData::gvalues)
        .def_readonly("genomes", &PyDiploidGeneticValueData::genomes);

    py::class_<genome_data_proxy>(m, "HaploidGenomeProxy")
        .def_readonly("effect_sizes", &genome_data_proxy::effect_sizes)
        .def_readonly("dominance", &genome_data_proxy::dominance)
        .def_readonly("positions", &genome_data_proxy::positions)
        .def_readonly("smutations", &genome_data_proxy::smutations);

    py::class_<double_array_proxy>(m, "_FloatProxy", py::buffer_protocol())
        .def_buffer([](const double_array_proxy& self) {
            return py::buffer_info(self.data, sizeof(double),
                                   py::format_descriptor<double>::format(), 1,
                                   {self.size}, {sizeof(double)});
        });

    py::class_<mutation_key_array_proxy>(m, "_MutationKeyArrayProxy",
                                         py::buffer_protocol())
        .def_buffer([](const mutation_key_array_proxy& self) {
            return py::buffer_info(self.data, sizeof(std::uint32_t),
                                   py::format_descriptor<std::uint32_t>::format(), 1,
                                   {self.size}, {sizeof(std::uint32_t)});
        });
}
