#include <boost/test/unit_test.hpp>

#include <core/evolve_discrete_demes/evolvets.hpp>
#include "fwdpy11/discrete_demography/DiscreteDemography.hpp"
#include "fwdpy11/discrete_demography/MigrationMatrix.hpp"
#include "fwdpy11/evolvets/SampleRecorder.hpp"
#include "fwdpy11/genetic_values/dgvalue_pointer_vector.hpp"
#include "fwdpy11/regions/MutationRegions.hpp"
#include <fwdpy11/genetic_values/fwdpp_wrappers/fwdpp_genetic_value.hpp>
#include <fwdpy11/genetic_value_to_fitness/GeneticValueIsTrait.hpp>
#include "fwdpy11/regions/RecombinationRegions.hpp"
#include "fwdpy11/rng.hpp"
#include "fwdpy11/types/DiploidPopulation.hpp"

BOOST_AUTO_TEST_SUITE(test_evolvets)

namespace
{
    // NOTE: all of this is a copy/paste
    // from code internal to the pybind11 module.

    struct single_deme_additive_het
    {
        inline void
        operator()(double& d, const fwdpy11::Mutation& m) const
        {
            d += m.s * m.h;
        }
    };

    struct multi_deme_additive_het
    {
        inline void
        operator()(const std::size_t deme, double& d, const fwdpy11::Mutation& m) const
        {
            d += m.esizes[deme] * m.heffects[deme];
        }
    };

    struct single_deme_additive_hom
    {
        double scaling;
        single_deme_additive_hom(double s) : scaling(s)
        {
        }

        inline void
        operator()(double& d, const fwdpy11::Mutation& m) const
        {
            d += scaling * m.s;
        }
    };

    struct multi_deme_additive_hom
    {
        double scaling;
        multi_deme_additive_hom(double s) : scaling(s)
        {
        }

        inline void
        operator()(const std::size_t deme, double& d, const fwdpy11::Mutation& m) const
        {
            d += scaling * m.esizes[deme];
        }
    };

    struct final_additive_trait
    {
        inline double
        operator()(double d) const
        {
            return d;
        }
    };

    struct final_additive_fitness
    {
        inline double
        operator()(double d) const
        {
            return std::max(0.0, 1. + d);
        }
    };

    using DiploidAdditive = fwdpy11::stateless_site_dependent_genetic_value_wrapper<
        single_deme_additive_het, single_deme_additive_hom, multi_deme_additive_het,
        multi_deme_additive_hom, 0>;
}

BOOST_AUTO_TEST_CASE(test_basic_api_coherence)
{
    fwdpy11::GSLrng_t rng(42);
    fwdpy11::DiploidPopulation pop(1000, 10.0);
    // no mutation
    fwdpy11::MutationRegions mregions({}, {});
    // no recombination
    fwdpy11::RecombinationRegions recregions(0., {});
    // no demographic events
    fwdpy11::discrete_demography::DiscreteDemography demography(
        {}, {}, {}, {}, fwdpy11::discrete_demography::MigrationMatrix(), {});
    DiploidAdditive additive(1, 2.0, final_additive_fitness(), nullptr, nullptr);

    fwdpy11::dgvalue_pointer_vector_ gvalue_ptrs(additive);

    // now, the callbacks...
    auto post_simplification_recorder = [](const fwdpy11::DiploidPopulation&) {};
    //auto stopping_criterion
    std::function<bool(const fwdpy11::DiploidPopulation&, const bool)> stopping_criterion
        = [](const fwdpy11::DiploidPopulation&, const bool) -> bool { return false; };
    auto sample_recorder_callback
        = [](const fwdpy11::DiploidPopulation&, const fwdpy11::SampleRecorder&) {};
    // Note: recorder cannot just be made internally
    // when the function is called because the memory
    // address must be known to pybind11.
    // In order to be passed to a Python-based recorder,
    // it must be created within Python.
    fwdpy11::SampleRecorder recorder;
    evolve_with_tree_sequences_options options;

    evolve_with_tree_sequences_refactor(rng, pop, recorder, 10, demography, 100, 0., 0.,
                                        mregions, recregions, gvalue_ptrs,
                                        sample_recorder_callback, stopping_criterion,
                                        post_simplification_recorder, options);
    BOOST_REQUIRE_EQUAL(pop.generation, 100);
}

BOOST_AUTO_TEST_SUITE_END()
