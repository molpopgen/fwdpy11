#ifndef FWDPY11_GSSMO
#define FWDPY11_GSSMO

#include <pybind11/stl.h>
#include <algorithm>
#include <vector>
#include <tuple>
#include "GeneticValueIsTrait.hpp"

namespace fwdpy11
{
    struct GSSmo : public GeneticValueIsTrait
    {
        double VS, opt;
        std::size_t current_optimum;
        // Tuple is time, optimum, VS
        std::vector<std::tuple<std::uint32_t, double, double>> optima;

        GSSmo(std::vector<std::tuple<std::uint32_t, double, double>> optima_)
            : GeneticValueIsTrait{ 1 },
              VS{ std::numeric_limits<double>::quiet_NaN() },
              opt{ std::numeric_limits<double>::quiet_NaN() },
              current_optimum(1), optima(std::move(optima_))
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

        double
        operator()(
            const DiploidMetadata &metadata,
            const std::vector<double> & /*genetic_values*/) const override
        {
            return std::exp(
                -(std::pow(metadata.g + metadata.e - opt, 2.0) / (2.0 * VS)));
        }

        template <typename poptype>
        inline void
        update_details(const poptype &pop)
        {
            if (current_optimum < optima.size())
                {
                    if (pop.generation >= std::get<0>(optima[current_optimum]))
                        {
                            opt = std::get<1>(optima[current_optimum]);
                            VS = std::get<2>(optima[current_optimum++]);
                        }
                }
        }

        void
        update(const DiploidPopulation &pop) override
        {
            update_details(pop);
        }

        std::unique_ptr<GeneticValueToFitnessMap>
        clone() const override
        {
            return std::unique_ptr<GSSmo>(new GSSmo(*this));
        }

        pybind11::object
        pickle() const override
        {
            return pybind11::make_tuple(opt, VS, current_optimum, optima);
        }
    };
} // namespace fwdpy11

#endif
