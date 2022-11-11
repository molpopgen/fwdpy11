#include <algorithm> //for std::max
#include <pybind11/pybind11.h>
#include <fwdpy11/genetic_values/DiploidGeneticValue.hpp>
#include <fwdpy11/genetic_values/default_update.hpp>

struct additive : public fwdpy11::DiploidGeneticValue
{
    additive() : fwdpy11::DiploidGeneticValue{1, nullptr, nullptr}
    {
    }

    double
    calculate_gvalue(const fwdpy11::DiploidGeneticValueData data) override
    {
        const auto& pop = data.pop.get();
        const auto diploid_index = data.offspring_metadata.get().label;
        double sum = 0;
        for (auto m : pop.haploid_genomes[pop.diploids[diploid_index].first].smutations)
            {
                sum += pop.mutations[m].s;
            }
        for (auto m : pop.haploid_genomes[pop.diploids[diploid_index].second].smutations)
            {
                sum += pop.mutations[m].s;
            }
        gvalues[0] = std::max(0.0, 1.0 + sum);
        return gvalues[0];
    }

    pybind11::object
    pickle() const
    {
        return pybind11::bytes("custom_additive");
    }
    DEFAULT_DIPLOID_POP_UPDATE();
};

//Standard pybind11 stuff goes here
PYBIND11_MODULE(custom_additive, m)
{
    pybind11::object imported_custom_additive_base_class_type
        = pybind11::module::import("fwdpy11").attr("DiploidGeneticValue");
    pybind11::class_<additive, fwdpy11::DiploidGeneticValue>(m, "additive")
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
