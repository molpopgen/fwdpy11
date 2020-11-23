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

import numpy as np
import pytest

import fwdpy11


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
    params = fwdpy11.ModelParams(**pdict)

    r = fwdpy11.RandomAncientSamples(seed=42, samplesize=2, timepoints=[3])

    fwdpy11.evolvets(rng, pop, params, 100, r)

    ts = pop.dump_tables_to_tskit()

    assert len(ts.tables.individuals) == pop.N + 2
    first = 0
    preserved = 0
    alive = 0
    for i in ts.tables.individuals:
        if i.flags & fwdpy11.tskit_tools.INDIVIDUAL_IS_ALIVE:
            alive += 1
            for n in i.metadata["nodes"]:
                assert ts.tables.nodes.time[n] == 0
        if i.flags & fwdpy11.tskit_tools.INDIVIDUAL_IS_PRESERVED:
            preserved += 1
            for n in i.metadata["nodes"]:
                assert (
                    ts.tables.nodes.time[n] == 0 + params.simlen
                    or ts.tables.nodes.time[n] == 0 + params.simlen - 3
                )

    assert alive == pop.N
    assert first == 0
    assert preserved == 2


@pytest.mark.parametrize("rng", [{"seed": 1234}], indirect=["rng"])
@pytest.mark.parametrize("pdict", [{"simlen": 10}], indirect=["pdict"])
@pytest.mark.parametrize("pop", [{"N": 100, "L": 1}], indirect=["pop"])
def test_metadata_roundtrip_single_sim_with_first_gen_preserved(rng, pdict, pop):
    params = fwdpy11.ModelParams(**pdict)

    r = fwdpy11.RandomAncientSamples(seed=42, samplesize=2, timepoints=[3])

    fwdpy11.evolvets(rng, pop, params, 100, r, preserve_first_generation=True)

    ts = pop.dump_tables_to_tskit()
    assert len(ts.tables.individuals) == 2 * pop.N + 2
    first = 0
    preserved = 0
    alive = 0
    for i in ts.tables.individuals:
        if i.flags & fwdpy11.tskit_tools.INDIVIDUAL_IS_ALIVE:
            alive += 1
            for n in i.metadata["nodes"]:
                assert ts.tables.nodes.time[n] == 0
        if i.flags & fwdpy11.tskit_tools.INDIVIDUAL_IS_PRESERVED:
            preserved += 1
            for n in i.metadata["nodes"]:
                assert (
                    ts.tables.nodes.time[n] == 0 + params.simlen
                    or ts.tables.nodes.time[n] == 0 + params.simlen - 3
                )
        if i.flags & fwdpy11.tskit_tools.INDIVIDUAL_IS_FIRST_GENERATION:
            first += 1
            for n in i.metadata["nodes"]:
                assert ts.tables.nodes.time[n] == 0 + params.simlen

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
        {"params_dict": str(params.asdict()), "script": inception}
    )

    provenance = json.loads(ts.provenance(0).record)
    params_dict = eval(provenance["parameters"]["params_dict"])

    assert params_dict == params.asdict()
    script = provenance["parameters"]["script"]

    with open(__file__, "r") as f:
        assert f.read() == script
