#include "pick_parents.hpp"
#include "fwdpy11/discrete_demography/exceptions.hpp"
#include <cwchar>
#include <iterator>
#include <stdexcept>

namespace fwdpy11_core
{
    namespace discrete_demography
    {
        parent_data
        pick_parents(const fwdpy11::GSLrng_t& rng, const std::int32_t offspring_deme,
                     const fwdpy11_core::ForwardDemesGraph& demography,
                     const fwdpp::gsl_ran_discrete_t_ptr& ancestor_deme_lookup,
                     const multideme_fitness_bookmark& fitness_bookmark,
                     const multideme_fitness_lookups<std::uint32_t>& wlookups,
                     bool allow_residual_selfing)
        {
            std::int32_t pdeme = static_cast<std::int32_t>(
                gsl_ran_discrete(rng.get(), ancestor_deme_lookup.get()));

            if (allow_residual_selfing == false)
                {
                    auto p = demography.parental_deme_sizes();
                    if (std::begin(p) == std::end(p))
                        {
                            throw std::runtime_error("parental deme sizes are NULL");
                        }
                    auto size = *(std::begin(p) + pdeme);
                    if (size == 1.0)
                        {
                            std::ostringstream o;
                            o << "residual selfing not allowed, but deme " << pdeme
                              << " has a size of 1";
                            throw fwdpy11::discrete_demography::DemographyError(o.str());
                        }
                }

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
            if (allow_residual_selfing == false && p1 == p2)
                {
                    while (p1 == p2)
                        {
                            p2 = wlookups.get_parent(rng, fitness_bookmark, pdeme);
                        }
                }
            return {p1, p2, pdeme, pdeme, mating_event_type::outcrossing};
        }
    }
}
