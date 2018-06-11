#ifndef FWDPY11_GENETIC_VALUE_TO_FITNESS_HPP__
#define FWDPY11_GENETIC_VALUE_TO_FITNESS_HPP__

#include <cmath>
#include <memory>
#include <functional>
#include <algorithm>
#include <limits>
#include <vector>
#include <queue>
#include <tuple>
#include <fwdpy11/types/SlocusPop.hpp>
#include <fwdpy11/types/MlocusPop.hpp>
#include <fwdpy11/genetic_values/default_update.hpp>

namespace fwdpy11
{
    using genetic_value_to_fitness_t
        = std::function<double(const double, const double)>;

    struct GeneticValueToFitness
    {
        virtual double operator()(const double, const double) const = 0;
        virtual void update(const SlocusPop &) = 0;
        virtual void update(const MlocusPop &) = 0;
        virtual std::unique_ptr<GeneticValueToFitness> clone() const = 0;
    };

    struct GeneticValueIsFitness : public GeneticValueToFitness
    {
        inline double
        operator()(const double g, const double) const
        {
            return g;
        }

        DEFAULT_SLOCUSPOP_UPDATE()
        DEFAULT_MLOCUSPOP_UPDATE()

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

        DEFAULT_SLOCUSPOP_UPDATE()
        DEFAULT_MLOCUSPOP_UPDATE()

        inline std::unique_ptr<GeneticValueToFitness>
        clone() const
        {
            return std::unique_ptr<GSS>(new GSS(VS, opt));
        }
    };

    struct GSSmo : public GeneticValueToFitness
    {
        double VS, opt;
        std::size_t current_optimum;
        // Tuple is time, optimum, VS
        std::vector<std::tuple<std::uint32_t, double, double>> optima;

        GSSmo(std::vector<std::tuple<std::uint32_t, double, double>>
                  optima_)
            : VS{ std::numeric_limits<double>::quiet_NaN() },
              opt{ std::numeric_limits<double>::quiet_NaN() }, current_optimum(1),optima(std::move(optima_))
        {
            using tuple_t = std::tuple<std::uint32_t, double, double>;
            if (optima.empty())
                {
                    throw std::invalid_argument("empty container of optima");
                }
            if (!std::is_sorted(optima.begin(), optima.end(),
                                [](const tuple_t &a, const tuple_t &b) {
                                    return std::get<0>(a) < std::get<0>(b);
                                }))
                {
                    throw std::invalid_argument("optima not sorted by time");
                }
            if (std::any_of(optima.begin(), optima.end(),
                            [](const tuple_t &t) {
                                auto VS_ = std::get<2>(t);
                                auto opt_ = std::get<1>(t);
                                bool rv = false;
                                if (VS_ < 0.0)
                                    rv = true;
                                if (!std::isfinite(VS_))
                                    rv = true;
                                if (!std::isfinite(opt_))
                                    rv = true;
                                return rv;
                            }))
                {
                    throw std::invalid_argument(
                        "all VS and opt values must be finite");
                }
            opt = std::get<1>(optima[0]);
            VS = std::get<2>(optima[0]);
        }

        inline double
        operator()(const double g, const double e) const
        {
            return std::exp(-(std::pow(g + e - opt, 2.0) / (2.0 * VS)));
        }

        template <typename poptype>
        inline void
        update_details(const poptype &pop)
        {
            if (current_optimum < optima.size())
                {
                    if (pop.generation >= std::get<0>(optima.front()))
                        {
                            opt = std::get<1>(optima[current_optimum]);
                            VS = std::get<2>(optima[current_optimum++]);
                        }
                }
        }

        inline void
        update(const SlocusPop &pop)
        {
            update_details(pop);
        }

        inline void
        update(const MlocusPop &pop)
        {
            update_details(pop);
        }

        inline std::unique_ptr<GeneticValueToFitness>
        clone() const
        {
            return std::unique_ptr<GSSmo>(new GSSmo(*this));
        }
    };
} //namespace fwdpy11

#endif
