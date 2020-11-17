#include <sstream>
#include <fwdpp/gsl_discrete.hpp>
#include <fwdpy11/discrete_demography/exceptions.hpp>
#include <fwdpy11/discrete_demography/simulation.hpp>
#include "discrete_demography_roundtrips.hpp"

MatingEventRecord::MatingEventRecord(std::uint32_t g, std::int32_t pd, std::int32_t od,
                                     fwdpy11::discrete_demography::mating_event_type mt)
    : generation{g}, parental_deme{pd}, offspring_deme{od}, mating_event{mt}
{
}

MockMetadata::MockMetadata(std::size_t parent1_, std::size_t parent2_,
                           std::size_t label_, std::int32_t deme_, double w_)
    : parent1{parent1_}, parent2{parent2_}, label{label_}, deme{deme_}, w{w_}
{
}

MockPopulation::MockPopulation(std::uint32_t popsize)

    : N{popsize}, generation{0}, diploid_metadata{init_metadata(popsize)}
{
}

std::vector<MockMetadata>
MockPopulation::init_metadata(std::uint32_t popsize)
{
    std::vector<MockMetadata> rv;
    for (decltype(popsize) i = 0; i < popsize; ++i)
        {
            rv.emplace_back(std::numeric_limits<std::size_t>::max(),
                            std::numeric_limits<std::size_t>::max(), static_cast<int>(i),
                            0, 1.0);
        }
    return rv;
}

std::unordered_map<int, int>
get_deme_sizes(const std::vector<MockMetadata>& metadata)
{
    std::unordered_map<int, int> rv;
    for (auto&& md : metadata)
        {
            rv[md.deme]++;
        }
    return rv;
}

std::vector<MatingEventRecord>
DiscreteDemography_roundtrip(
    const fwdpy11::GSLrng_t& rng, MockPopulation& pop,
    fwdpy11::discrete_demography::DiscreteDemography& demography, std::uint32_t ngens)
// Assumptions:
// 1. Initial deme labels are set by the user. NOTE: validated by manager object
//
{
    auto current_demographic_state
        = fwdpy11::discrete_demography::initialize_model_state(
            pop.generation, pop.diploid_metadata, demography);
    decltype(pop.diploid_metadata) offspring_metadata;
    offspring_metadata.reserve(pop.N);
    std::vector<MatingEventRecord> rv;
    fwdpy11::discrete_demography::update_demography_manager(
        rng, pop.generation, pop.diploid_metadata, demography,
        *current_demographic_state);
    if (current_demographic_state->will_go_globally_extinct() == true)
        {
            std::ostringstream o;
            o << "extinction at time " << pop.generation;
            throw fwdpy11::discrete_demography::GlobalExtinction(o.str());
        }
    for (std::uint32_t gen = 0; gen < ngens; ++gen)
        {
            ++pop.generation;
            // Generate the offspring
            for (std::int32_t deme = 0; deme < current_demographic_state->maxdemes;
                 ++deme)
                {
                    for (unsigned ind = 0; ind < current_demographic_state->sizes_rates
                                                     .next_deme_sizes.get()[deme];
                         ++ind)
                        {
                            auto pdata = fwdpy11::discrete_demography::pick_parents(
                                rng, deme, current_demographic_state->miglookup,
                                current_demographic_state->sizes_rates
                                    .current_deme_sizes,
                                current_demographic_state->sizes_rates.selfing_rates,
                                current_demographic_state->fitnesses);
                            if (pdata.deme1 != deme)
                                {
                                    rv.emplace_back(pop.generation + 1, pdata.deme1,
                                                    deme, pdata.mating);
                                }
                            if (pdata.deme2 != deme)
                                {
                                    rv.emplace_back(pop.generation + 1, pdata.deme2,
                                                    deme, pdata.mating);
                                }
                            offspring_metadata.emplace_back(pdata.parent1, pdata.parent2,
                                                            offspring_metadata.size(),
                                                            deme, 1.);
                        }
                }
            pop.diploid_metadata.swap(offspring_metadata);
            offspring_metadata.clear();

            fwdpy11::discrete_demography::update_demography_manager(
                rng, pop.generation, pop.diploid_metadata, demography,
                *current_demographic_state);
            if (current_demographic_state->will_go_globally_extinct() == true)
                {
                    std::ostringstream o;
                    o << "extinction at time " << pop.generation;
                    throw fwdpy11::discrete_demography::GlobalExtinction(o.str());
                }
            // NOTE: this test doesn't populate the diploid genotypes.
            // Only the metadata "evolve", so this update differs
            // from what we do in a "real" simulation.
            pop.N = static_cast<std::uint32_t>(pop.diploid_metadata.size());
        }
    fwdpy11::discrete_demography::save_model_state(std::move(current_demographic_state),
                                                   demography);
    return rv;
}

