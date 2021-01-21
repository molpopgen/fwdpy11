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

import numpy as np
import pytest

import fwdpy11


class CountSamplesPerTimePoint(object):
    def __init__(self):
        self.sample_timepoints = []
        self.sample_sizes = []
        self.timepoint_seen = {}

    def __call__(self, pop):
        assert len(pop.tables.preserved_nodes) // 2 == len(pop.ancient_sample_metadata)
        # Get the most recent ancient samples
        # and record their number.  We do this
        # by a "brute-force" approach
        for t, n, m in pop.sample_timepoints(False):
            if t not in self.timepoint_seen:
                self.timepoint_seen[t] = 1
            else:
                self.timepoint_seen[t] += 1
            if t not in self.sample_timepoints:
                self.sample_timepoints.append(t)
                self.sample_sizes.append(len(n) // 2)

            # simplify to each time point
            tables, idmap = fwdpy11.simplify_tables(pop.tables, n)
            for ni in n:
                assert idmap[ni] != fwdpy11.NULL_NODE
                assert tables.nodes[idmap[ni]].time == t


@pytest.fixture
def params(request):
    pd = {
        "nregions": request.param["nregions"],
        "sregions": [fwdpy11.ExpS(beg=0, end=1, weight=1, mean=0.2)],
        "recregions": [fwdpy11.PoissonInterval(0, 1, 1e-2)],
        "gvalue": [fwdpy11.Multiplicative(2.0)],
        "rates": request.param["rates"],
        "demography": fwdpy11.DiscreteDemography(),
        "simlen": request.param["simlen"],
    }
    return fwdpy11.ModelParams(**pd)


@pytest.fixture
def pop():
    return fwdpy11.DiploidPopulation(1000, 1.0)


@pytest.mark.parametrize("seed", [101 * 45 * 110 * 210])
@pytest.mark.parametrize(
    "params",
    [
        {"simlen": 101, "nregions": [], "rates": (0.0, 0.025, None)},
        {
            "simlen": 101,
            "nregions": [fwdpy11.Region(0.0, 1.0, 1.0)],
            "rates": (0.025, 0.025, None),
        },
    ],
    indirect=["params"],
)
def test_ancient_sample_resetting(pop, params, seed):
    rng = fwdpy11.GSLrng(seed)
    ancient_sample_recorder = fwdpy11.RandomAncientSamples(
        seed=42, samplesize=10, timepoints=[i for i in range(1, 101)]
    )
    resetter = CountSamplesPerTimePoint()
    fwdpy11.evolvets(
        rng,
        pop,
        params,
        5,
        recorder=ancient_sample_recorder,
        post_simplification_recorder=resetter,
    )
    assert len(pop.tables.preserved_nodes) == 0
    assert len(pop.ancient_sample_metadata) == 0

    for _, j in resetter.timepoint_seen.items():
        assert j == 1

    assert resetter.sample_timepoints == [i for i in range(1, 101)]

    assert all([i == 10 for i in resetter.sample_sizes]) is True
