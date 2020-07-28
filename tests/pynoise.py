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
import fwdpy11
import fwdpy11.custom_genetic_value_decorators


@fwdpy11.custom_genetic_value_decorators.default_update
@fwdpy11.custom_genetic_value_decorators.genetic_value_noise_default_clone
class PyNoise(fwdpy11.GeneticValueNoise):
    def __init__(self):
        fwdpy11.GeneticValueNoise.__init__(self)

    def __call__(self, data: fwdpy11.DiploidGeneticValueNoiseData) -> float:
        return fwdpy11.gsl_ran_gaussian_ziggurat(data.rng, 0.1)
