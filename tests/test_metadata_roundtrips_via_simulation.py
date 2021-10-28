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

import json

import fwdpy11
import msprime
import numpy as np
import pytest


@pytest.fixture
def pdict(request):
    pd = {
        "recregions": [],
        "sregions": [],
        "nregions": [],
        "rates": (0, 0, 0),
        "gvalue": fwdpy11.Multiplicative(2.0),
        "demography": fwdpy11.DiscreteDemography(),
        "simlen": request.param["simlen"],
    }
    return pd


@pytest.fixture
def pop(request):
    return fwdpy11.DiploidPopulation(request.param["N"], request.param["L"])


@pytest.fixture
def inception():
    with open(__file__, "r") as f:
        contents = f.read()
    return contents


@pytest.mark.parametrize("rng", [{"seed": 666}], indirect=["rng"])
@pytest.mark.parametrize("pdict", [{"simlen": 10}], indirect=["pdict"])
@pytest.mark.parametrize("pop", [{"N": 100, "L": 1}], indirect=["pop"])
def test_metadata_roundtrip_single_sim(rng, pdict, pop):
    import sys

    params = fwdpy11.ModelParams(**pdict)

    r = fwdpy11.RandomAncientSamples(seed=42, samplesize=2, timepoints=[3])

    fwdpy11.evolvets(rng, pop, params, 100, r)

    ts = pop.dump_tables_to_tskit(model_params=params)

    # add neutral mutations w/no metadata
    ts = msprime.sim_mutations(ts, rate=1.0, random_seed=654321)
    # bulk decode the mutation metadata, which is all None
    mutation_md = fwdpy11.tskit_tools.decode_mutation_metadata(ts)
    assert all([i is None for i in mutation_md])

    # test index access
    num_muts = 0
    for i in range(ts.num_mutations):
        x = fwdpy11.tskit_tools.decode_mutation_metadata(ts, i)
        assert x[0] is None
        num_muts += 1
    assert num_muts == ts.num_mutations

    num_muts = 0
    for i in fwdpy11.tskit_tools.decode_mutation_metadata(
        ts, slice(0, ts.num_mutations, 2)
    ):
        assert i is None
        num_muts += 1

    assert num_muts == len([i for i in range(0, ts.num_mutations, 2)])

    assert ts.num_individuals == pop.N + 2
    first = 0
    preserved = 0
    alive = 0
    for i in ts.individuals():
        if i.flags & fwdpy11.tskit_tools.INDIVIDUAL_IS_ALIVE:
            alive += 1
            for n in i.metadata["nodes"]:
                assert ts.node(n).time == 0
        if i.flags & fwdpy11.tskit_tools.INDIVIDUAL_IS_PRESERVED:
            preserved += 1
            for n in i.metadata["nodes"]:
                assert (
                    ts.node(n).time == 0 + params.simlen
                    or ts.node(n).time == 0 + params.simlen - 3
                )

    assert alive == pop.N
    assert first == 0
    assert preserved == 2
    # Test slice/index access to individual metadata
    first = 0
    preserved = 0
    alive = 0
    for row in range(ts.num_individuals):
        i = fwdpy11.tskit_tools.decode_individual_metadata(ts, row)[0]
        assert type(i.alive) == bool
        assert type(i.preserved) == bool
        assert type(i.first_generation) == bool
        if i.alive:
            alive += 1
            for n in i.nodes:
                assert ts.node(n).time == 0
        if i.preserved:
            preserved += 1
            for n in i.nodes:
                assert (
                    ts.node(n).time == 0 + params.simlen
                    or ts.node(n).time == 0 + params.simlen - 3
                )
    assert alive == pop.N
    assert first == 0
    assert preserved == 2

    num_md = 0
    for _ in fwdpy11.tskit_tools.decode_individual_metadata(
        ts, slice(1, ts.num_individuals, 2)
    ):
        num_md += 1
    assert num_md == len([i for i in range(0, ts.num_individuals, 2)])


@pytest.mark.parametrize("rng", [{"seed": 1234}], indirect=["rng"])
@pytest.mark.parametrize("pdict", [{"simlen": 10}], indirect=["pdict"])
@pytest.mark.parametrize("pop", [{"N": 100, "L": 1}], indirect=["pop"])
def test_metadata_roundtrip_single_sim_with_first_gen_preserved(rng, pdict, pop):
    params = fwdpy11.ModelParams(**pdict)

    r = fwdpy11.RandomAncientSamples(seed=42, samplesize=2, timepoints=[3])

    fwdpy11.evolvets(rng, pop, params, 100, r, preserve_first_generation=True)

    ts = pop.dump_tables_to_tskit(model_params=params)
    assert ts.num_individuals == 2 * pop.N + 2
    first = 0
    preserved = 0
    alive = 0
    for i in ts.individuals():
        if i.flags & fwdpy11.tskit_tools.INDIVIDUAL_IS_ALIVE:
            alive += 1
            for n in i.metadata["nodes"]:
                assert ts.node(n).time == 0
        if i.flags & fwdpy11.tskit_tools.INDIVIDUAL_IS_PRESERVED:
            preserved += 1
            for n in i.metadata["nodes"]:
                assert (
                    ts.node(n).time == 0 + params.simlen
                    or ts.node(n).time == 0 + params.simlen - 3
                )
        if i.flags & fwdpy11.tskit_tools.INDIVIDUAL_IS_FIRST_GENERATION:
            first += 1
            for n in i.metadata["nodes"]:
                assert ts.node(n).time == 0 + params.simlen

    assert alive == pop.N
    assert first == pop.N
    assert preserved == 2 + pop.N


@pytest.mark.parametrize("rng", [{"seed": 666}], indirect=["rng"])
@pytest.mark.parametrize("pdict", [{"simlen": 10}], indirect=["pdict"])
@pytest.mark.parametrize("pop", [{"N": 100, "L": 1}], indirect=["pop"])
def test_metadata_roundtrip_single_deme_sim_with_parameters(rng, pdict, pop, inception):
    params = fwdpy11.ModelParams(**pdict)

    fwdpy11.evolvets(rng, pop, params, 100)

    ts = pop.dump_tables_to_tskit(
        parameters={"params_dict": str(params.asdict()), "script": inception}
    )

    provenance = json.loads(ts.provenance(0).record)
    params_dict = eval(provenance["parameters"]["params_dict"])

    assert params_dict == params.asdict()
    script = provenance["parameters"]["script"]

    with open(__file__, "r") as f:
        assert f.read() == script
