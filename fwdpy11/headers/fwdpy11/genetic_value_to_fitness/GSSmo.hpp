#ifndef FWDPY11_GSSMO
#define FWDPY11_GSSMO

#include <algorithm>
#include <vector>
#include "GeneticValueIsTrait.hpp"
#include "Optimum.hpp"

namespace fwdpy11
{
    struct GSSmo : public GeneticValueIsTrait
    {
        double VS, opt;
        std::size_t current_optimum;
        // Tuple is time, optimum, VS
        std::vector<Optimum> optima;

        GSSmo(std::vector<Optimum> optima_)
            : GeneticValueIsTrait{1}, VS{std::numeric_limits<double>::quiet_NaN()},
              opt{std::numeric_limits<double>::quiet_NaN()}, current_optimum(1),
              optima(std::move(optima_))
        {
            if (optima.empty())
                {
                    throw std::invalid_argument("empty container of optima");
                }
            for (auto &o : optima)
                {
                    if (o.when == Optimum::null)
                        {
                            throw std::invalid_argument(
                                "invalid when value for Optimum");
                        }
                }
            if (!std::is_sorted(
                    optima.begin(), optima.end(),
                    [](const Optimum &a, const Optimum &b) { return a.when < b.when; }))
                {
                    throw std::invalid_argument("optima not sorted by time");
                }
            opt = optima[0].opt;
            VS = optima[0].VW;
        }

        double
        operator()(const DiploidGeneticValueToFitnessData data) const override
        {
            return std::exp(-(std::pow(data.offspring_metadata.get().g
                                           + data.offspring_metadata.get().e - opt,
                                       2.0)
                              / (2.0 * VS)));
        }

        template <typename poptype>
        inline void
        update_details(const poptype &pop)
        {
            if (current_optimum < optima.size())
                {
                    if (pop.generation >= optima[current_optimum].when)
                        {
                            opt = optima[current_optimum].opt;
                            VS = optima[current_optimum++].VW;
                        }
                }
        }

        void
        update(const DiploidPopulation &pop) override
        {
            update_details(pop);
        }

        std::shared_ptr<GeneticValueToFitnessMap>
        clone() const override
        {
            return std::make_shared<GSSmo>(optima);
        }
    };
} // namespace fwdpy11

#endif
