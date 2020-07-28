#
# Copyright (C) 2017-2020 Kevin Thornton <krthornt@uci.edu>
#
# This file is part of fwdpy11.
#
# fwdpy11 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# fwdpy11 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with fwdpy11.  If not, see <http://www.gnu.org/licenses/>.
#

import math

import fwdpy11
import fwdpy11.custom_genetic_value_decorators


@fwdpy11.custom_genetic_value_decorators.default_update
class PyAdditiveGSS(fwdpy11.PyDiploidGeneticValue):
    def __init__(self, opt, VS):
        self.opt = opt
        self.VS = VS
        fwdpy11.PyDiploidGeneticValue.__init__(self, 1, None, None)

    def genetic_value_to_fitness(
        self, data: fwdpy11.DiploidGeneticValueToFitnessData
    ) -> float:
        return math.e ** (
            -((data.offspring_metadata.g + data.offspring_metadata.e - self.opt) ** 2)
            / (2 * self.VS)
        )

    def calculate_gvalue(self, data: fwdpy11.PyDiploidGeneticValueData) -> float:
        s = fwdpy11.strict_additive_effects(data.pop, data.offspring_metadata)
        memoryview(data)[0] = s
        return s
