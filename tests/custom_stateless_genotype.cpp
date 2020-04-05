#include <algorithm> //for std::max
#include <pybind11/pybind11.h>
#include <fwdpp/fitness_models.hpp>
#include <fwdpy11/genetic_values/DiploidGeneticValue.hpp>
#include <fwdpy11/genetic_values/default_update.hpp>

struct GeneralW : public fwdpy11::DiploidGeneticValue
{
    fwdpp::site_dependent_genetic_value w;

    GeneralW() : fwdpy11::DiploidGeneticValue{ 1 }, w{} {}

    inline double
    calculate_gvalue(const std::size_t diploid_index,
                     const fwdpy11::DiploidMetadata& /*metadata*/,
                     const fwdpy11::DiploidPopulation& pop) const override
    {
        gvalues[0]
            = std::max(0.0, w(
                                pop.diploids[diploid_index],
                                pop.haploid_genomes, pop.mutations,
                                [](double& g, const fwdpy11::Mutation& m) {
                                    g *= (1.0 + m.s);
                                },
                                [](double& g, const fwdpy11::Mutation& m) {
                                    g *= (1.0 + m.h);
                                }, 1.0));
        return gvalues[0];
    }

    pybind11::object
    pickle() const override
    {
        return pybind11::bytes("custom_stateless_genotype");
    }

    DEFAULT_DIPLOID_POP_UPDATE();
};

//Standard pybind11 stuff goes here
PYBIND11_MODULE(custom_stateless_genotype, m)
{
    pybind11::object imported_custom_stateless_genotype_base_class_type
        = pybind11::module::import("fwdpy11").attr("DiploidGeneticValue");
    pybind11::class_<GeneralW, fwdpy11::DiploidGeneticValue>(m, "GeneralW")
        .def(pybind11::init<>())
        .def(pybind11::pickle([](const GeneralW& g) { return g.pickle(); },
                              [](pybind11::object o) {
                                  pybind11::bytes b(o);
                                  auto s = b.cast<std::string>();
                                  if (s != "custom_stateless_genotype")
                                      {
                                          throw std::runtime_error(
                                              "invalid object state");
                                      }
                                  return GeneralW();
                              }));
}
