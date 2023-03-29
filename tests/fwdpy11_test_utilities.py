#
# Copyright (C) 2021 Kevin Thornton <krthornt@uci.edu>
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


def seed_list(seed, nseeds):
    rng = np.random.Generator(np.random.MT19937(seed=seed))
    used = {}
    seeds = []
    for _ in range(nseeds):
        candidate = rng.integers(0, np.iinfo(np.uint32).max)
        while candidate in used:
            candidate = rng.integers(0, np.iinfo(np.uint32).max)
        used[candidate] = 1
        seeds.append(candidate)
    return seeds
