#include "core/demes/forward_graph.hpp"
#include "functions.hpp"
#include "fwdpp/internal/sample_diploid_helpers.hpp"
#include <cstddef>
#include <cstdint>

namespace
{
    void
    parent_deme_exists(
        fwdpy11_core::ForwardDemesGraphDataIterator<double> parental_deme_sizes,
        std::uint32_t parental_deme, std::size_t offspring_deme,
        std::uint32_t generation)
    {
        auto parental_deme_size = *(std::begin(parental_deme_sizes) + parental_deme);
        if (parental_deme_size == 0.0)
            {
                std::ostringstream o;
                o << "deme " << offspring_deme
                  << " has ancestry from "
                     "deme "
                  << parental_deme << " at time " << generation
                  << " but the parental "
                     "deme size "
                     "is 0";
                throw fwdpy11::discrete_demography::DemographyError(o.str());
            }
    }

    void
    fitnesses_exist(const fwdpy11_core::discrete_demography::multideme_fitness_lookups<
                        std::uint32_t>& fitnesses,
                    std::uint32_t parental_deme, std::uint32_t generation)
    {
        if (fitnesses.lookups[parental_deme] == nullptr)
            {
                std::ostringstream o;
                o << "fitness lookup table "
                     "for "
                     "parental deme "
                  << parental_deme << " is empty at time " << generation;
                throw fwdpy11::discrete_demography::DemographyError(o.str());
            }
    }

    void
    process_deme(std::uint32_t generation, std::size_t offspring_deme,
                 fwdpy11_core::ForwardDemesGraphDataIterator<double> parental_deme_sizes,
                 const fwdpy11_core::discrete_demography::multideme_fitness_lookups<
                     std::uint32_t>& fitnesses,
                 const fwdpy11_core::ForwardDemesGraph& demography)
    { // offspring deme exists
        auto ancestry_proportions
            = demography.offspring_ancestry_proportions(offspring_deme);
        std::uint32_t parental_deme = 0;
        for (auto prop = std::begin(ancestry_proportions);
             prop != std::end(ancestry_proportions); ++prop, ++parental_deme)
            {
                if (*prop > 0.0)
                    {
                        parent_deme_exists(parental_deme_sizes, parental_deme,
                                           offspring_deme, generation);
                        fitnesses_exist(fitnesses, parental_deme, generation);
                    }
            }
    }

}

namespace fwdpy11_core
{
    namespace discrete_demography
    {
        void
        validate_parental_state(
            std::uint32_t generation,
            const multideme_fitness_lookups<std::uint32_t>& fitnesses,
            const fwdpy11_core::ForwardDemesGraph& demography)
        {
            if (demography.iterating_model())
                {
                    auto offspring_deme_sizes = demography.offspring_deme_sizes();
                    auto parental_deme_sizes = demography.parental_deme_sizes();

                    std::size_t offspring_deme = 0;
                    for (auto size = std::begin(offspring_deme_sizes);
                         size != std::end(offspring_deme_sizes);
                         ++size, ++offspring_deme)
                        {
                            if (*size > 0.0)
                                {
                                    process_deme(generation, offspring_deme,
                                                 parental_deme_sizes, fitnesses,
                                                 demography);
                                }
                        }
                }
        }
    }
}
