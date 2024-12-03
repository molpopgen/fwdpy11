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

import demes
import fwdpy11
import pytest
import tskit

from dataclasses import dataclass


@pytest.fixture
def pdict1():
    pd = {
        "nregions": [],
        "sregions": [fwdpy11.ExpS(0, 1, 1, -0.2)],
        "recregions": [fwdpy11.PoissonInterval(0, 1, 1e-3)],
        "gvalue": [fwdpy11.Additive(2.0)],
        "rates": (0, 1e-4, None),
        "simlen": 1000,
        "prune_selected": False,
    }
    return pd


@pytest.fixture
def pdict2():
    pd = {
        "nregions": [],
        "sregions": [fwdpy11.ExpS(0, 1, 1, -0.2)],
        "recregions": [],
        "gvalue": [fwdpy11.Additive(2.0)],
        "rates": (0, 1e-4, None),
        "simlen": 1000,
        "prune_selected": False,
    }
    return pd


@pytest.mark.parametrize("pop", [{"N": 100, "genome_length": 1}], indirect=["pop"])
def test_single_model_params(pop, pdict1):
    demography = fwdpy11.ForwardDemesGraph.tubes([100], 1)
    pdict1["demography"] = demography
    mp = fwdpy11.ModelParams(**pdict1)

    ts = pop.dump_tables_to_tskit(model_params=mp)

    # reconstruct
    for i in dir(fwdpy11):
        exec(f"from fwdpy11 import {i}")

    mp_rebuilt = fwdpy11.ModelParams(**eval(ts.metadata["model_params"]))

    assert mp == mp_rebuilt


@pytest.mark.parametrize("pop", [{"N": 100, "genome_length": 1}], indirect=["pop"])
def test_population_metadata(pop, gutenkunst):
    dm = fwdpy11.ForwardDemesGraph.from_demes(gutenkunst, burnin=10)
    demes_metadata = {}
    for key, value in dm.deme_labels.items():
        demes_metadata[key] = {"name": value}
    ts = pop.dump_tables_to_tskit(population_metadata=demes_metadata)

    for i in range(ts.tables.populations.num_rows):
        assert ts.population(i).metadata == demes_metadata[i]

    # Method 2: get the description out, too
    demes_metadata = {}
    for i, d in enumerate(gutenkunst.demes):
        demes_metadata[i] = {"name": d.name, "description": d.description}

    # Validate that dm thinks that things are in the same order
    for key, value in dm.deme_labels.items():
        assert demes_metadata[key]["name"] == value

    ts = pop.dump_tables_to_tskit(population_metadata=demes_metadata)
    for i in range(ts.tables.populations.num_rows):
        md = ts.population(i).metadata
        assert md["name"] == gutenkunst.demes[i].name
        assert md["description"] == gutenkunst.demes[i].description


@pytest.mark.parametrize("pop", [{"N": 100, "genome_length": 1}], indirect=["pop"])
def test_user_defined_data(pop):
    ts = pop.dump_tables_to_tskit(data={"mydata": 11})
    assert ts.metadata["data"]["mydata"] == 11

    class MyType(object):
        def __init__(self, x):
            self.x = x

        def __repr__(self):
            return f"MyType(x={self.x})"

    ts = pop.dump_tables_to_tskit(data=str(MyType(x=11)))
    assert eval(ts.metadata["data"]).x == 11


@pytest.mark.parametrize("pop", [{"N": 100, "genome_length": 1}], indirect=["pop"])
def test_seed(pop):
    ts = pop.dump_tables_to_tskit(seed=333)
    assert ts.metadata["seed"] == 333

    with pytest.raises(ValueError):
        _ = pop.dump_tables_to_tskit(seed=-333)


def test_toplevel_metadata_0():
    tc = tskit.TableCollection(10.0)
    tc.metadata_schema = fwdpy11.tskit_tools.metadata_schema.TopLevelMetadata
    md = {"generation": 100, "data": {0.0: [0, 1, 2]}}
    tc.metadata = md
    data = tc.metadata["data"]
    assert 0.0 not in data
    assert "0.0" in data


def test_toplevel_metadata_1():
    # Same as previous test, but w/no schema
    # and we just shove raw bytes in there
    tc = tskit.TableCollection(10.0)
    md = {"generation": 100, "data": {0.0: [0, 1, 2]}}
    tc.metadata = str(md).encode()
    decoded = eval(tc.metadata)
    data = decoded["data"]
    assert 0.0 in data


def test_toplevel_metadata_2():
    # This time w/a JSON codec
    tc = tskit.TableCollection(10.0)
    schema = tskit.metadata.MetadataSchema(
        {
            "codec": "json",
        }
    )
    tc.metadata_schema = schema
    md = {"generation": 100, "data": {0.0: [0, 1, 2]}}

    tc.metadata = md
    decoded = tc.metadata
    data = decoded["data"]
    assert 0.0 not in data
    assert "0.0" in data


@dataclass(eq=True, frozen=True)
class FloatWrap:
    value: float


def test_toplevel_metadata_3():

    # This time w/a JSON codec
    tc = tskit.TableCollection(10.0)
    schema = tskit.metadata.MetadataSchema(
        {
            "codec": "json",
        }
    )
    tc.metadata_schema = schema
    md = {"data": {FloatWrap(0.0): [0, 1, 2]}}
    tc.metadata = str(md)
    decoded = eval(tc.metadata)
    data = decoded["data"]
    assert FloatWrap(0.0) in data
