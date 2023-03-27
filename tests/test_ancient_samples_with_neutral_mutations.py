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
from hypothesis import settings, given, HealthCheck
from hypothesis.strategies import booleans, integers

import pytest

import fwdpy11


class CountSamplesPerTimePoint(object):
    def __init__(self):
        self.sample_timepoints = []
        self.sample_sizes = []
        self.timepoint_seen = {}

    def __call__(self, pop):
        # Get the most recent ancient samples
        # and record their number.  We do this
        # by a "brute-force" approach
        for t, n, _ in pop.sample_timepoints(False):
            assert t != pop.generation
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


@pytest.fixture(
    scope="function",
    params=[
        {"nregions": [], "rates": (0.0, 0.025, None)},
        {
            "nregions": [fwdpy11.Region(0.0, 1.0, 1.0)],
            "rates": (0.025, 0.025, None),
        },
    ],
)
def pdict(request):
    """
    Build models w/ and w/o neutral mutations.
    """
    pd = {
        "nregions": request.param["nregions"],
        "sregions": [fwdpy11.ExpS(beg=0, end=1, weight=1, mean=0.2)],
        "recregions": [fwdpy11.PoissonInterval(0, 1, 1e-2)],
        "gvalue": [fwdpy11.Multiplicative(2.0)],
        "rates": request.param["rates"],
        "demography": None,
    }
    return pd


def make_counter():
    return lambda: CountSamplesPerTimePoint()


def make_none():
    return lambda: None


@pytest.fixture(scope="function", params=[make_counter, make_none])
def resetter(request):
    """
    Generate a value for post_simplification_recorder
    """
    return request.param()()


@pytest.fixture(params=[5, 10, 25, 33, 50, 73, 100, 101])
def simlen(request):
    return request.param


@pytest.fixture(params=[1, 2, 5, 10, 33])
def simplification_inteval(request):
    return request.param


@pytest.fixture(params=[True, False])
def prune_selected(request):
    return request.param


@given(popsize=integers(50, 200),
       simlen=integers(5, 102),
       simplification_inteval=integers(1, 50),
       seed=integers(1, int(2**32 - 1)),
       prune_selected=booleans(),
       resetter=booleans())
@settings(suppress_health_check=[HealthCheck.function_scoped_fixture],
          deadline=None, max_examples=5)
def test_ancient_samples_and_neutral_mutations(
    pdict, simlen, simplification_inteval, prune_selected, seed, popsize, resetter
):
    """
    The test involving neutral mutations test GitHub issue 639 and 646
    """
    mslike_pop = fwdpy11.DiploidPopulation(popsize, 1.0)
    demography = fwdpy11.ForwardDemesGraph.tubes(mslike_pop.deme_sizes()[1],
                                                 burnin=simlen,
                                                 burnin_is_exact=True)
    if resetter is False:
        ancient_sample_resetter = None
    else:
        ancient_sample_resetter = CountSamplesPerTimePoint()
    rng = fwdpy11.GSLrng(seed)
    ancient_sample_recorder = fwdpy11.RandomAncientSamples(
        seed=42, samplesize=10, timepoints=[i for i in range(1, simlen + 1)]
    )
    pdict["simlen"] = demography.final_generation
    pdict["prune_selected"] = prune_selected
    pdict["demography"] = demography
    params = fwdpy11.ModelParams(**pdict)
    fwdpy11.evolvets(
        rng,
        mslike_pop,
        params,
        simplification_inteval,
        recorder=ancient_sample_recorder,
        post_simplification_recorder=ancient_sample_resetter,
    )

    if ancient_sample_resetter is not None:
        assert len(mslike_pop.ancient_sample_metadata) == 0

        # FIXME: check that the number of sampled time points is correct?
        for _, j in ancient_sample_resetter.timepoint_seen.items():
            assert j == 1

        assert all(
            [i == 10 for i in ancient_sample_resetter.sample_sizes]) is True
    else:
        assert (
            len(mslike_pop.ancient_sample_metadata)
            == len([i for i in range(1, params.simlen)]) * 10
        )

    ti = fwdpy11.TreeIterator(
        mslike_pop.tables, mslike_pop.alive_nodes, update_samples=True
    )
    for t in ti:
        for m in t.mutations():
            s = t.samples_below(m.node)
            if len(s) == 0:
                assert mslike_pop.mcounts_ancient_samples[m.key] > 0
            else:
                assert len(s) == mslike_pop.mcounts[m.key]
                if prune_selected is True and len(s) == 2 * mslike_pop.N:
                    assert mslike_pop.mcounts_ancient_samples[m.key] > 0
