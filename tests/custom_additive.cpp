#include <algorithm> //for std::max
#include <pybind11/pybind11.h>
#include <fwdpy11/genetic_values/DiploidPopulationGeneticValue.hpp>
#include <fwdpy11/genetic_values/default_update.hpp>

    struct additive : public fwdpy11::DiploidPopulationGeneticValue
{
    additive() : fwdpy11::DiploidPopulationGeneticValue(1) {}

    inline double
    calculate_gvalue(const std::size_t diploid_index,
                     const fwdpy11::DiploidPopulation& pop) const
    {
        double sum = 0;
        for (auto m :
             pop.haploid_genomes[pop.diploids[diploid_index].first].smutations)
            {
                sum += pop.mutations[m].s;
            }
        for (auto m :
             pop.haploid_genomes[pop.diploids[diploid_index].second].smutations)
            {
                sum += pop.mutations[m].s;
            }
        gvalues[0] = std::max(0.0, 1.0 + sum);
        return gvalues[0];
    }

    double
    genetic_value_to_fitness(const fwdpy11::DiploidMetadata& metadata) const
    {
        return metadata.g;
    }

    double
    noise(const fwdpy11::GSLrng_t& /*rng*/,
          const fwdpy11::DiploidMetadata& /*offspring_metadata*/,
          const std::size_t /*parent1*/, const std::size_t /*parent2*/,
          const fwdpy11::DiploidPopulation& /*pop*/) const
    {
        return 0.0;
    }

    pybind11::object
    pickle() const
    {
        return pybind11::bytes("custom_additive");
    }
    DEFAULT_DIPLOID_POP_UPDATE();

    pybind11::tuple
    shape() const
    {
        return pybind11::make_tuple(1);
    }
};

//Standard pybind11 stuff goes here
PYBIND11_MODULE(custom_additive, m)
{
    pybind11::object imported_custom_additive_base_class_type
        = pybind11::module::import("fwdpy11")
              .attr("GeneticValue");
    pybind11::class_<additive, fwdpy11::DiploidPopulationGeneticValue>(m, "additive")
        .def(pybind11::init<>())
        .def(pybind11::pickle([](const additive& a) { return a.pickle(); },
                              [](pybind11::object o) {
                                  pybind11::bytes b(o);
                                  auto s = b.cast<std::string>();
                                  if (s != "custom_additive")
                                      {
                                          throw std::runtime_error(
                                              "invalid object state");
                                      }
                                  return additive();
                              }));
}
