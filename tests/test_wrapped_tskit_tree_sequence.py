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

import fwdpy11
import msprime
import numpy as np
import pytest


@pytest.fixture(scope="session")
def simulation():
    pdict = {
        "nregions": [],
        "sregions": [fwdpy11.ExpS(0, 1, 1, -1e-3)],
        "recregions": [fwdpy11.PoissonInterval(0, 1, 1e-3)],
        "rates": (0, 1e-3, None),
        "gvalue": fwdpy11.Multiplicative(2.0),
        "simlen": 100,
    }
    pop = fwdpy11.DiploidPopulation(1000, 1.0)
    rng = fwdpy11.GSLrng(101)
    recorder = fwdpy11.RandomAncientSamples(101, 10, np.arange(1, pdict["simlen"]))

    params = fwdpy11.ModelParams(**pdict)

    fwdpy11.evolvets(rng, pop, params, 100, recorder=recorder)

    return pop, params


@pytest.mark.parametrize("pop", [{"N": 100, "genome_length": 1}], indirect=["pop"])
def test_creation(pop):
    ts = pop.dump_tables_to_tskit()

    ts.dump("foo.trees")

    fwdpy11.tskit_tools.WrappedTreeSequence(ts)


def test_on_simulated_example(simulation):
    pop, params = simulation
    ts = pop.dump_tables_to_tskit(model_params=params)

    wts = fwdpy11.tskit_tools.WrappedTreeSequence(ts=ts)

    assert wts.model_params == params

    times = []
    for time, nodes, md in wts.timepoints_with_individuals(decode_metadata=True):
        assert np.all(wts.ts.tables.nodes.time[nodes] == time)
        for i in md:
            if time == 0.0:
                assert i.alive is True
                assert i.preserved is False
            else:
                assert i.alive is False
                assert i.preserved is True
            for n in i.nodes:
                assert wts.ts.tables.nodes.time[n] == time
        times.append(time)

    assert np.all(np.array(times) == np.arange(pop.generation, dtype=np.float64)[::-1])

    ## Test storing two model param object
    ts = pop.dump_tables_to_tskit(model_params={"stage1": params, "stage2": params})
    wts = fwdpy11.tskit_tools.WrappedTreeSequence(ts)
    mp = wts.model_params
    assert len(mp) == 2
    for _, value in mp.items():
        assert value == params


def test_generation_property(simulation):
    pop, _ = simulation
    assert pop.generation > 0
    ts = pop.dump_tables_to_tskit()
    wts = fwdpy11.tskit_tools.WrappedTreeSequence(ts=ts)
    assert wts.generation == pop.generation


def test_demes_graph_property(gutenkunst):
    pop = fwdpy11.DiploidPopulation(100, 1.0)
    ts = pop.dump_tables_to_tskit(demes_graph=gutenkunst)
    wts = fwdpy11.tskit_tools.WrappedTreeSequence(ts)
    assert wts.demes_graph == gutenkunst
    ts = pop.dump_tables_to_tskit()
    wts = fwdpy11.tskit_tools.WrappedTreeSequence(ts)
    assert wts.demes_graph is None


def test_with_msprime_provenance_added():
    initial_ts = msprime.sim_ancestry(samples=200, population_size=100, random_seed=1)
    pop = fwdpy11.DiploidPopulation.create_from_tskit(initial_ts)
    ts = pop.dump_tables_to_tskit()
    ts = msprime.sim_mutations(ts, rate=1e-5, random_seed=1)
    assert ts.tables.provenances.num_rows == 2
    wts = fwdpy11.tskit_tools.WrappedTreeSequence(ts)
    assert wts.ts.tables.provenances.num_rows == 2
