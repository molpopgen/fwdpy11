#include <boost/test/unit_test.hpp>

#include <core/evolve_discrete_demes/evolvets.hpp>
#include <core/demes/forward_graph.hpp>
#include <fwdpy11/evolvets/SampleRecorder.hpp>
#include <fwdpy11/genetic_values/dgvalue_pointer_vector.hpp>
#include <fwdpy11/regions/MutationRegions.hpp>
#include <fwdpy11/genetic_values/fwdpp_wrappers/fwdpp_genetic_value.hpp>
#include <fwdpy11/genetic_value_to_fitness/GeneticValueIsTrait.hpp>
#include <fwdpy11/regions/RecombinationRegions.hpp>
#include <fwdpy11/rng.hpp>
#include <fwdpy11/types/DiploidPopulation.hpp>

#include "forward_demes_graph_fixtures.hpp"

BOOST_AUTO_TEST_SUITE(test_evolvets)

/* Tests needed
 *
 * - [X] simlen > model time in graph
 * - [] initial population config not compatible
 *      with state of parental demes
 * - [] simlen < model time in graph,
 *      but we keep simulating until we are done.
 * - [] One deme, multiple epochs, test size hitory is 
 *      correct
 */

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

struct common_setup
{
    fwdpy11::GSLrng_t rng;
    fwdpy11::DiploidPopulation pop;
    fwdpy11::MutationRegions mregions;
    fwdpy11::RecombinationRegions recregions;
    DiploidAdditive additive;
    fwdpy11::dgvalue_pointer_vector_ gvalue_ptrs;
    fwdpy11::SampleRecorder recorder;
    std::function<void(const fwdpy11::DiploidPopulation&)> post_simplification_recorder;
    std::function<void(const fwdpy11::DiploidPopulation&,
                       const fwdpy11::SampleRecorder&)>
        sample_recorder_callback;
    std::function<bool(const fwdpy11::DiploidPopulation&, const bool)>
        stopping_criterion;
    evolve_with_tree_sequences_options options;

    common_setup()
        : rng{42}, pop{100, 10.0}, mregions{{}, {}}, recregions{0., {}},
          additive{1, 2.0, final_additive_fitness(), nullptr, nullptr},
          gvalue_ptrs{additive}, recorder{},
          post_simplification_recorder{[](const fwdpy11::DiploidPopulation&) {}},
          sample_recorder_callback{
              [](const fwdpy11::DiploidPopulation&, const fwdpy11::SampleRecorder&) {}},
          stopping_criterion{[](const fwdpy11::DiploidPopulation&, const bool) -> bool {
              return false;
          }},
          options{}
    {
    }
};

BOOST_FIXTURE_TEST_CASE(test_basic_api_coherence, common_setup)
{
    auto model = SingleDemeModel();
    fwdpy11_core::ForwardDemesGraph forward_demes_graph(model.yaml, 10);

    // TODO: if we put long run times in here, we get exceptions
    // from the ForwardDemesGraph back end.
    evolve_with_tree_sequences_refactor(rng, pop, recorder, 10, forward_demes_graph, 10,
                                        0., 0., mregions, recregions, gvalue_ptrs,
                                        sample_recorder_callback, stopping_criterion,
                                        post_simplification_recorder, options);
    BOOST_REQUIRE_EQUAL(pop.generation, 10);
}

BOOST_FIXTURE_TEST_CASE(test_simlen_longer_than_model_length, common_setup)
{
    auto model = SingleDemeModel();
    fwdpy11_core::ForwardDemesGraph forward_demes_graph(model.yaml, 10);

    // Here, simlen = 100, but the graph only has 10 generations of births in it.
    evolve_with_tree_sequences_refactor(rng, pop, recorder, 10, forward_demes_graph, 100,
                                        0., 0., mregions, recregions, gvalue_ptrs,
                                        sample_recorder_callback, stopping_criterion,
                                        post_simplification_recorder, options);
    BOOST_REQUIRE_EQUAL(pop.generation, 10);
}

BOOST_AUTO_TEST_SUITE_END()
