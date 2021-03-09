
//
// Copyright (C) 2017-2020 Kevin Thornton <krthornt@uci.edu>
//
// This file is part of fwdpy11.
//
// fwdpy11 is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// fwdpy11 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with fwdpy11.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef FWDPY11_GENETIC_MAP_UNIT_HPP
#define FWDPY11_GENETIC_MAP_UNIT_HPP

#include "fwdpp/genetic_map/fixed_number_crossovers.hpp"
#include <memory>
#include <cstdint>
#include <utility>
#include <fwdpp/genetic_map/genetic_map_unit.hpp>
#include <fwdpp/genetic_map/poisson_interval.hpp>
#include <fwdpp/genetic_map/poisson_point.hpp>
#include <fwdpp/genetic_map/binomial_interval.hpp>
#include <fwdpp/genetic_map/binomial_point.hpp>

namespace fwdpy11
{
    class GeneticMapUnit
    /// Base class To hold the low-level fwdpp
    /// types.  The python interface will create them,
    /// and then new instances will be cloned off to the C++ side.
    {
      private:
        std::unique_ptr<fwdpp::genetic_map_unit> ll_map_unit;

        bool discrete_;

      public:
        using ll_ptr_t = std::unique_ptr<fwdpp::genetic_map_unit>;
        explicit GeneticMapUnit(ll_ptr_t&& input, bool discrete)
            : ll_map_unit(std::move(input)), discrete_(discrete)
        {
        }
        GeneticMapUnit(GeneticMapUnit&&)=default;

        GeneticMapUnit(const GeneticMapUnit &) = delete;
        GeneticMapUnit& operator=(const GeneticMapUnit &) = delete;
        GeneticMapUnit& operator=(GeneticMapUnit &&) = default;


        ll_ptr_t
        ll_clone() const
        {
            return ll_map_unit->clone();
        }

        virtual ~GeneticMapUnit() = default;

        bool
        discrete() const
        {
            return discrete_;
        }
    };

    template <typename fwdpp_type_double, typename fwdpp_type_discrete>
    struct GeneticMapUnitInitializationTrampoline
    /// Trampoline class manages the boiler plate
    {
      public:
        GeneticMapUnitInitializationTrampoline()
        {
        }

        template <typename... Args>
        inline GeneticMapUnit::ll_ptr_t
        operator()(bool discrete, Args&&... args) const
        {
            if (discrete)
                {
                    return std::make_unique<fwdpp_type_discrete>(
                        std::forward<Args>(args)...);
                }
            return std::make_unique<fwdpp_type_double>(std::forward<Args>(args)...);
        }
    };

    // Derived classes hold pointers to either continuous (double)
    // or integer-valued (int64_t) objects.

    struct PoissonInterval : GeneticMapUnit
    {
        using trampoline_t = GeneticMapUnitInitializationTrampoline<
            fwdpp::poisson_interval, fwdpp::poisson_interval_t<std::int64_t>>;
        PoissonInterval(double begin, double end, double mean, bool discrete)
            : GeneticMapUnit(trampoline_t()(discrete, begin, end, mean), discrete)
        {
        }
    };

    struct PoissonPoint : public GeneticMapUnit
    {
      public:
        using trampoline_t = GeneticMapUnitInitializationTrampoline<
            fwdpp::poisson_point, fwdpp::poisson_point_t<std::int64_t>>;
        PoissonPoint(double position, double mean, bool discrete)
            : GeneticMapUnit(trampoline_t()(discrete, position, mean), discrete)
        {
        }
    };

    struct BinomialInterval : GeneticMapUnit
    {
        using trampoline_t = GeneticMapUnitInitializationTrampoline<
            fwdpp::binomial_interval, fwdpp::binomial_interval_t<std::int64_t>>;
        BinomialInterval(double begin, double end, double probability, bool discrete)
            : GeneticMapUnit(trampoline_t()(discrete, begin, end, probability), discrete)
        {
        }
    };

    struct BinomialPoint : public GeneticMapUnit
    {
      public:
        using trampoline_t = GeneticMapUnitInitializationTrampoline<
            fwdpp::binomial_point, fwdpp::binomial_point_t<std::int64_t>>;
        BinomialPoint(double position, double probability, bool discrete)
            : GeneticMapUnit(trampoline_t()(discrete, position, probability), discrete)
        {
        }
    };

    struct FixedCrossovers : public GeneticMapUnit
    {
        using trampoline_t = GeneticMapUnitInitializationTrampoline<
            fwdpp::fixed_number_crossovers,
            fwdpp::fixed_number_crossovers_t<std::int64_t>>;
        FixedCrossovers(double begin, double end, int nxovers, bool discrete)
            : GeneticMapUnit(trampoline_t()(discrete, begin, end, nxovers), discrete)
        {
        }
    };
}

#endif
