#include <functional>
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <fwdpp/fitness_models.hpp>
#include <fwdpy11/genetic_values/GeneticValueToFitness.hpp>
#include <fwdpy11/genetic_values/SlocusPopGeneticValue.hpp>

namespace py = pybind11;

struct single_locus_additive : public fwdpy11::SlocusPopGeneticValue
{
    const fwdpp::additive_diploid gv;
    const fwdpy11::genetic_value_to_fitness_t gv2w;

    single_locus_additive(const double scaling)
        : gv{ scaling, fwdpp::additive_diploid::policy::aw }, gv2w{
            fwdpy11::GeneticValueIsFitness()
          }
    {
    }

    single_locus_additive(const double scaling,
                          fwdpp::additive_diploid::policy p,
                          const fwdpy11::genetic_value_to_fitness_t& g2w)
        : gv{ scaling, p }, gv2w{ g2w }
    {
    }

    inline double
    operator()(const std::size_t diploid_index,
               const fwdpy11::SlocusPop& pop) const
    {
        return gv(pop.diploids[diploid_index], pop.gametes, pop.mutations);
    }

    inline void
    update(const fwdpy11::SlocusPop&)
    {
    }

    inline double
    genetic_value_to_fitness(const double g, const double e) const
    {
        return gv2w(g, e);
    }
};

PYBIND11_MODULE(genetic_values, m)
{
    py::enum_<fwdpp::additive_diploid::policy>(m, "AdditivePolicy")
        .value("aw", fwdpp::additive_diploid::policy::aw)
        .value("atrait", fwdpp::additive_diploid::policy::atrait)
        .export_values();

    py::class_<single_locus_additive>(m, "SlocusAdditive")
        .def(py::init<double>(), py::arg("scaling"))
        .def(py::init<double, fwdpp::additive_diploid::policy,
                      fwdpy11::genetic_value_to_fitness_t>());
}
