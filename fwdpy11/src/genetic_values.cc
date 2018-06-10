#include <functional>
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <fwdpp/fitness_models.hpp>
#include <fwdpy11/genetic_values/GeneticValueToFitness.hpp>
#include <fwdpy11/genetic_values/SlocusPopGeneticValue.hpp>

namespace py = pybind11;

template <typename fwdppT>
struct wrap_fwdpp_genetic_value : public fwdpy11::SlocusPopGeneticValue
{
    const fwdppT gv;
    const fwdpy11::genetic_value_to_fitness_t gv2w;

    wrap_fwdpp_genetic_value(const double);

    wrap_fwdpp_genetic_value(const double scaling,
                             const typename fwdppT::policy p,
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

template <>
wrap_fwdpp_genetic_value<fwdpp::additive_diploid>::wrap_fwdpp_genetic_value(
    const double scaling)
    : gv{ scaling, fwdpp::additive_diploid::policy::aw }, gv2w{
          fwdpy11::GeneticValueIsFitness()
      }
{
}

template <>
wrap_fwdpp_genetic_value<fwdpp::multiplicative_diploid>::
    wrap_fwdpp_genetic_value(const double scaling)
    : gv{ scaling, fwdpp::multiplicative_diploid::policy::mw }, gv2w{
          fwdpy11::GeneticValueIsFitness()
      }
{
}

using wrapped_additive = wrap_fwdpp_genetic_value<fwdpp::additive_diploid>;
using wrapped_multiplicative
    = wrap_fwdpp_genetic_value<fwdpp::multiplicative_diploid>;

PYBIND11_MODULE(genetic_values, m)
{
    py::enum_<fwdpp::additive_diploid::policy>(m, "AdditivePolicy")
        .value("aw", fwdpp::additive_diploid::policy::aw)
        .value("atrait", fwdpp::additive_diploid::policy::atrait)
        .export_values();

    py::enum_<fwdpp::multiplicative_diploid::policy>(m, "MultiplicativePolicy")
        .value("mw", fwdpp::multiplicative_diploid::policy::mw)
        .value("mtrait", fwdpp::multiplicative_diploid::policy::mtrait)
        .export_values();

    py::class_<wrapped_additive>(m, "SlocusAdditive")
        .def(py::init<double>(), py::arg("scaling"))
        .def(py::init<double, fwdpp::additive_diploid::policy,
                      fwdpy11::genetic_value_to_fitness_t>())
        .def_property_readonly(
            "scaling",
            [](const wrapped_additive& wa) { return wa.gv.scaling; })
        .def_property_readonly("is_fitness", [](const wrapped_additive& wa) {
            return wa.gv.p == fwdpp::additive_diploid::policy::aw;
        });

    py::class_<wrapped_multiplicative>(m, "SlocusMult")
        .def(py::init<double>(), py::arg("scaling"))
        .def(py::init<double, fwdpp::multiplicative_diploid::policy,
                      fwdpy11::genetic_value_to_fitness_t>())
        .def_property_readonly(
            "scaling",
            [](const wrapped_multiplicative& wa) { return wa.gv.scaling; })
        .def_property_readonly(
            "is_fitness", [](const wrapped_multiplicative& wa) {
                return wa.gv.p == fwdpp::multiplicative_diploid::policy::mw;
            });
}
