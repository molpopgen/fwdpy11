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

    return pop


@pytest.mark.parametrize("pop", [{"N": 100, "genome_length": 1}], indirect=["pop"])
def test_creation(pop):
    ts = pop.dump_tables_to_tskit()

    ts.dump("foo.trees")

    fwdpy11.tskit_tools.WrappedTreeSequence(ts)


def test_on_simulated_example(simulation):
    ts = simulation.dump_tables_to_tskit()

    wts = fwdpy11.tskit_tools.WrappedTreeSequence(ts=ts)

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

    assert np.all(
        np.array(times) == np.arange(simulation.generation, dtype=np.float64)[::-1]
    )
