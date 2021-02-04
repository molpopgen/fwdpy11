#
# Copyright (C) 2017 Kevin Thornton <krthornt@uci.edu>
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
import pytest

import fwdpy11


@pytest.fixture(scope="function")
def rng(request):
    try:
        seed = request.param.get("seed", None)
        return fwdpy11.GSLrng(seed)
    except AttributeError:  # NOQA
        return fwdpy11.GSLrng(42)


@pytest.fixture(scope="function")
def pop(request):
    popsize = request.param["N"]
    genome_length = request.param["genome_length"]
    return fwdpy11.DiploidPopulation(popsize, genome_length)


@pytest.fixture(scope="function")
def mslike_pop(request):
    try:
        N = request.param.get("N", None)
        return fwdpy11.DiploidPopulation(N, 1.0)
    except AttributeError:  # NOQA
        return fwdpy11.DiploidPopulation(1000, 1.0)
