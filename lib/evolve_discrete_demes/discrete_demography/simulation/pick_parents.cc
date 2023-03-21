#include "pick_parents.hpp"

namespace fwdpy11_core
{
    namespace discrete_demography
    {
        parent_data
        pick_parents(const fwdpy11::GSLrng_t& rng, const std::int32_t offspring_deme,
                     const fwdpy11_core::ForwardDemesGraph& demography,
                     const fwdpp::gsl_ran_discrete_t_ptr& ancestor_deme_lookup,
                     const multideme_fitness_bookmark& fitness_bookmark,
                     const multideme_fitness_lookups<std::uint32_t>& wlookups)
        {
            std::int32_t pdeme = static_cast<std::int32_t>(
                gsl_ran_discrete(rng.get(), ancestor_deme_lookup.get()));

            auto selfing_rates = demography.offspring_selfing_rates();
            auto selfing_rate_offspring_deme
                = *(std::begin(selfing_rates) + offspring_deme);

            auto p1 = wlookups.get_parent(rng, fitness_bookmark, pdeme);
            // FIXME: this gives rise to residual selfing
            if (selfing_rate_offspring_deme > 0.
                && gsl_rng_uniform(rng.get()) <= selfing_rate_offspring_deme)
                {
                    return {p1, p1, pdeme, pdeme, mating_event_type::selfing};
                }
            auto p2 = wlookups.get_parent(rng, fitness_bookmark, pdeme);
            return {p1, p2, pdeme, pdeme, mating_event_type::outcrossing};
        }
    }
}
