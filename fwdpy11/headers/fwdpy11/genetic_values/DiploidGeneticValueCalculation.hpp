//
// Copyright (C) 2025 Kevin Thornton <krthornt@uci.edu>
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

#pragma once

#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpy11/genetic_value_data/genetic_value_data.hpp>

namespace fwdpy11
{
    struct DiploidGeneticValueCalculation
    {
        virtual ~DiploidGeneticValueCalculation() = default;

        virtual double calculate_gvalue(const DiploidGeneticValueData data) = 0;
        virtual void update(const DiploidPopulation& pop) = 0;
        virtual std::size_t ndim() const = 0;
        virtual std::shared_ptr<DiploidGeneticValueCalculation> clone() const = 0;
    };
}
