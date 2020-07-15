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

import numpy as np

import fwdpy11


class PyAdditive(fwdpy11.PyDiploidGeneticValue):
    def __init__(self, gvalue_to_fitness=None, noise=None):
        fwdpy11.PyDiploidGeneticValue.__init__(self, 1, gvalue_to_fitness, noise)

    def calculate_gvalue(self, data):
        s = 0.0
        for g in data.genomes:
            s += np.array(g.effect_sizes, copy=False).sum()
        np.array(data.gvalues, copy=False)[0] = s
        return s
