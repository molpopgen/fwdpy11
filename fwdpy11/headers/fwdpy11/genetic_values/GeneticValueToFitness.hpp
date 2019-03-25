//
// Copyright (C) 2017 Kevin Thornton <krthornt@uci.edu>
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
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpy11/genetic_values/default_update.hpp>

namespace fwdpy11
{
    struct GeneticValueToFitnessMap
    {
        virtual ~GeneticValueToFitnessMap() = default;
        virtual double
        operator()(const DiploidMetadata & /*metadata*/) const = 0;
        virtual void update(const DiploidPopulation & /*pop*/) = 0;
        virtual std::unique_ptr<GeneticValueToFitnessMap> clone() const = 0;
        virtual pybind11::object pickle() const = 0;
    };

    struct GeneticValueIsFitness : public GeneticValueToFitnessMap
    {
        inline double
        operator()(const DiploidMetadata &metadata) const
        {
            return metadata.g;
        }

        DEFAULT_DIPLOID_POP_UPDATE()

        inline std::unique_ptr<GeneticValueToFitnessMap>
        clone() const
        {
            return std::unique_ptr<GeneticValueIsFitness>(
                new GeneticValueIsFitness());
        }

        virtual pybind11::object
        pickle() const
        {
            return pybind11::bytes("GeneticValueIsFitness");
        }
    };

    struct GeneticValueIsTrait : public GeneticValueToFitnessMap
    /// Another ABC.  Effectively a type trait
    {
    };

    struct GSS : public GeneticValueIsTrait
    {
        const double opt, VS;
        GSS(const double opt_, const double VS_) : opt{ opt_ }, VS{ VS_ }
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
        operator()(const DiploidMetadata &metadata) const
        {
            return std::exp(
                -(std::pow(metadata.g + metadata.e - opt, 2.0) / (2.0 * VS)));
        }

        DEFAULT_DIPLOID_POP_UPDATE()

        inline std::unique_ptr<GeneticValueToFitnessMap>
        clone() const
        {
            return std::unique_ptr<GSS>(new GSS(opt, VS));
        }

        virtual pybind11::object
        pickle() const
        {
            return pybind11::make_tuple(opt, VS);
        }
    };

    struct GSSmo : public GeneticValueIsTrait
    {
        double VS, opt;
        std::size_t current_optimum;
        // Tuple is time, optimum, VS
        std::vector<std::tuple<std::uint32_t, double, double>> optima;

        GSSmo(std::vector<std::tuple<std::uint32_t, double, double>> optima_)
            : VS{ std::numeric_limits<double>::quiet_NaN() },
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

        inline double
        operator()(const DiploidMetadata &metadata) const
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

        inline void
        update(const DiploidPopulation &pop)
        {
            update_details(pop);
        }

        inline std::unique_ptr<GeneticValueToFitnessMap>
        clone() const
        {
            return std::unique_ptr<GSSmo>(new GSSmo(*this));
        }

        virtual pybind11::object
        pickle() const
        {
            return pybind11::make_tuple(opt, VS, current_optimum, optima);
        }
    };
} //namespace fwdpy11

#endif
