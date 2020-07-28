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

import attr
import numpy as np

import fwdpy11
import fwdpy11.custom_genetic_value_decorators


@fwdpy11.custom_genetic_value_decorators.default_update
@fwdpy11.custom_genetic_value_decorators.genetic_value_is_trait_default_clone()
@attr.s()
class PyGSS(fwdpy11.GeneticValueIsTrait):
    opt = attr.ib()
    VS = attr.ib()

    def __attrs_post_init__(self):
        fwdpy11.GeneticValueIsTrait.__init__(self)

    def __call__(self, data: fwdpy11.DiploidGeneticValueToFitnessData) -> float:
        return math.e ** (
            -((data.offspring_metadata.g + data.offspring_metadata.e - self.opt) ** 2)
            / (2 * self.VS)
        )


@fwdpy11.custom_genetic_value_decorators.genetic_value_is_trait_default_clone()
@attr.s()
class PyGSSRandomOptimum(fwdpy11.GeneticValueIsTrait):
    opt = attr.ib()
    VS = attr.ib()

    def __attrs_post_init__(self):
        fwdpy11.GeneticValueIsTrait.__init__(self)
        self.optima = []

    def __call__(self, data: fwdpy11.DiploidGeneticValueToFitnessData) -> float:
        return math.e ** (
            -((data.offspring_metadata.g + data.offspring_metadata.e - self.opt) ** 2)
            / (2 * self.VS)
        )

    def update(self, pop: fwdpy11.DiploidPopulation) -> None:
        self.opt = np.random.normal(0.0, 1.0, 1,)[0]
        self.optima.append((pop.generation, self.opt))
