#ifndef FWDPY11_GENETIC_VALUE_TO_FITNESS_HPP__
#define FWDPY11_GENETIC_VALUE_TO_FITNESS_HPP__

#include <cmath>
#include <memory>
#include <functional>
#include <fwdpy11/types/SlocusPop.hpp>
#include <fwdpy11/types/MlocusPop.hpp>

namespace fwdpy11
{
    using genetic_value_to_fitness_t
        = std::function<double(const double, const double)>;

    struct GeneticValueToFitness
    {
        inline virtual double operator()(const double, const double) const = 0;
        inline virtual void update(const SlocusPop &) = 0;
        inline virtual void update(const MlocusPop &) = 0;
        inline virtual std::unique_ptr<GeneticValueToFitness>
        clone() const = 0;
    };

    struct GeneticValueIsFitness : public GeneticValueToFitness
    {
        inline double
        operator()(const double g, const double) const
        {
            return g;
        }
        inline void
        update(const SlocusPop &)
        {
        }
        inline void
        update(const MlocusPop &)
        {
        }
        inline std::unique_ptr<GeneticValueToFitness>
        clone() const
        {
            return std::unique_ptr<GeneticValueIsFitness>(
                new GeneticValueIsFitness());
        }
    };

    struct GSS : public GeneticValueToFitness
    {
        const double VS, opt;
        GSS(const double VS_, const double opt_) : VS{ VS_ }, opt{ opt_ }
        {
            if (VS <= 0.0)
                {
                    throw std::invalid_argument("VS must be > 0.0");
                }
            if (!std::isfinite(VS) || !std::isfinite(opt))
                {
                    throw std::invalid_argument(
                        "Both VS and opt must be finite values");
                }
        }

        inline double
        operator()(const double g, const double e) const
        {
            return std::exp(-(std::pow(g + e - opt, 2.0) / (2.0 * VS)));
        }

        inline void
        update(const SlocusPop &)
        {
        }

        inline void
        update(const MlocusPop &)
        {
        }

        inline std::unique_ptr<GeneticValueToFitness>
        clone() const
        {
            return std::unique_ptr<GSS>(new GSS(VS, opt));
        }
    };
} //namespace fwdpy11

#endif
