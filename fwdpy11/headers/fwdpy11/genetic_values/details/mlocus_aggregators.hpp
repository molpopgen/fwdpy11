#ifndef FWDPY11_GENETIC_VALUES_DETAILS_MLOCUS_AGGREGATORS_HPP__
#define FWDPY11_GENETIC_VALUES_DETAILS_MLOCUS_AGGREGATORS_HPP__

#include <vector>
#include <numeric>
#include <functional>

namespace fwdpy11
{
    struct aggregate_additive_fitness
    {
        inline double
        operator()(const std::vector<double>& g) const noexcept
        {
            auto s = g.size();
            return std::max(0., std::accumulate(g.data(), g.data() + s, 0.0)
                                    - (s - 1));
        }
    };

    struct aggregate_mult_fitness
    {
        inline double
        operator()(const std::vector<double>& g) const noexcept
        {
            auto s = g.size();
            return std::max(0., std::accumulate(g.data(), g.data() + s, 1.0,
                                                std::multiplies<double>()));
        }
    };

    struct aggregate_additive_trait
    {
        inline double
        operator()(const std::vector<double>& g) const noexcept
        {
            return std::accumulate(g.data(), g.data() + g.size(), 0.0);
        }
    };

    struct aggregate_mult_trait
    {
        inline double
        operator()(const std::vector<double>& g) const noexcept
        {
            auto s = g.size();
            return std::accumulate(
                       g.data(), g.data() + s, 1.0,
                       [](double prod, double v) { return prod * (1. + v); })
                   - 1.0;
        }
    };
} // namespace fwdpy11

#endif
