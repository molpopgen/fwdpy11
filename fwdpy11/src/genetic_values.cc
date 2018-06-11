#include <cmath>
#include <memory>
#include <functional>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <fwdpp/fitness_models.hpp>
#include <fwdpy11/genetic_values/GeneticValueToFitness.hpp>
#include <fwdpy11/genetic_values/SlocusPopGeneticValue.hpp>
#include <fwdpy11/genetic_values/default_update.hpp>

namespace py = pybind11;

template <typename fwdppT>
struct wrap_fwdpp_genetic_value : public fwdpy11::SlocusPopGeneticValue
{
    const fwdppT gv;
    const std::unique_ptr<fwdpy11::GeneticValueToFitness> gv2w;

    wrap_fwdpp_genetic_value(const double);

    wrap_fwdpp_genetic_value(const double scaling,
                             const typename fwdppT::policy p,
                             const fwdpy11::GeneticValueToFitness& g2w)
        : gv{ scaling, p }, gv2w{ g2w.clone() }
    {
        if (!std::isfinite(scaling))
            {
                throw std::invalid_argument("scaling must be finite");
            }
    }

    inline double
    operator()(const std::size_t diploid_index,
               const fwdpy11::SlocusPop& pop) const
    {
        return gv(pop.diploids[diploid_index], pop.gametes, pop.mutations);
    }

    DEFAULT_SLOCUSPOP_UPDATE()

    inline double
    genetic_value_to_fitness(const double g, const double e) const
    {
        return gv2w->operator()(g, e);
    }
};

template <>
wrap_fwdpp_genetic_value<fwdpp::additive_diploid>::wrap_fwdpp_genetic_value(
    const double scaling)
    : gv{ scaling, fwdpp::additive_diploid::policy::aw }, gv2w{
          std::unique_ptr<fwdpy11::GeneticValueIsFitness>(
              new fwdpy11::GeneticValueIsFitness())
      }
{
    if (!std::isfinite(scaling))
        {
            throw std::invalid_argument("scaling must be finite");
        }
}

template <>
wrap_fwdpp_genetic_value<fwdpp::multiplicative_diploid>::
    wrap_fwdpp_genetic_value(const double scaling)
    : gv{ scaling, fwdpp::multiplicative_diploid::policy::mw }, gv2w{
          std::unique_ptr<fwdpy11::GeneticValueIsFitness>(
              new fwdpy11::GeneticValueIsFitness())
      }
{
}

using wrapped_additive = wrap_fwdpp_genetic_value<fwdpp::additive_diploid>;
using wrapped_multiplicative
    = wrap_fwdpp_genetic_value<fwdpp::multiplicative_diploid>;

struct GBR : public fwdpy11::SlocusPopGeneticValue
{
    const std::unique_ptr<fwdpy11::GeneticValueToFitness> gv2w;
    GBR(const fwdpy11::GeneticValueToFitness& g2w) : gv2w{ g2w.clone() } {}

    inline double
    sum_haplotype_effect_sizes(const std::size_t gamete_index,
                               const fwdpy11::SlocusPop& pop) const
    {
        double h = 0.0;
        for (auto&& key : pop.gametes[gamete_index].smutations)
            {
                h += pop.mutations[key].s;
            }
        return h;
    }

    inline double
    operator()(const std::size_t diploid_index,
               const fwdpy11::SlocusPop& pop) const
    {
        double h1 = sum_haplotype_effect_sizes(
            pop.diploids[diploid_index].first, pop);
        double h2 = sum_haplotype_effect_sizes(
            pop.diploids[diploid_index].second, pop);
        return std::sqrt(h1 * h2);
    }

    DEFAULT_SLOCUSPOP_UPDATE()

    inline double
    genetic_value_to_fitness(const double g, const double e) const
    {
        return gv2w->operator()(g, e);
    }
};

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

    py::class_<fwdpy11::SlocusPopGeneticValue>(
        m, "SlocusPopGeneticValue",
        "ABC for genetic value calculations for diploid members of "
        ":class:`fwdpy11.SlocusPop`");

    py::class_<wrapped_additive, fwdpy11::SlocusPopGeneticValue>(
        m, "SlocusAdditive")
        .def(py::init<double>(), py::arg("scaling"))
        .def(py::init<double, fwdpp::additive_diploid::policy,
                      const fwdpy11::GeneticValueToFitness&>())
        .def_property_readonly(
            "scaling",
            [](const wrapped_additive& wa) { return wa.gv.scaling; })
        .def_property_readonly("is_fitness", [](const wrapped_additive& wa) {
            return wa.gv.p == fwdpp::additive_diploid::policy::aw;
        });

    py::class_<wrapped_multiplicative, fwdpy11::SlocusPopGeneticValue>(
        m, "SlocusMult")
        .def(py::init<double>(), py::arg("scaling"))
        .def(py::init<double, fwdpp::multiplicative_diploid::policy,
                      const fwdpy11::GeneticValueToFitness&>())
        .def_property_readonly(
            "scaling",
            [](const wrapped_multiplicative& wa) { return wa.gv.scaling; })
        .def_property_readonly(
            "is_fitness", [](const wrapped_multiplicative& wa) {
                return wa.gv.p == fwdpp::multiplicative_diploid::policy::mw;
            });

    py::class_<GBR, fwdpy11::SlocusPopGeneticValue>(m, "GBR").def(
        py::init<const fwdpy11::GeneticValueToFitness&>(), py::arg("gv2w"));

    py::class_<fwdpy11::GeneticValueToFitness>(
        m, "GeneticValueToFitness",
        "ABC for functions translating genetic values into fitness.");

    py::class_<fwdpy11::GeneticValueIsFitness, fwdpy11::GeneticValueToFitness>(
        m, "GeneticValueIsFitness")
        .def(py::init<>());

    py::class_<fwdpy11::GSS, fwdpy11::GeneticValueToFitness>(
        m, "GSS", "Gaussian stabilizing selection.")
        .def(py::init<double, double>(), py::arg("VS"), py::arg("opt"));

    py::class_<fwdpy11::GSSmo, fwdpy11::GeneticValueToFitness>(
        m, "GSSmo", "Gaussian stabilizing selection with a moving optimum.")
        .def(
            py::init<std::vector<std::tuple<std::uint32_t, double, double>>>(),
            py::arg("optima"));
}
