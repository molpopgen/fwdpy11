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

import warnings

import demes
import fwdpy11
import pytest


@pytest.fixture
def pdict1():
    pd = {
        "nregions": [],
        "sregions": [fwdpy11.ExpS(0, 1, 1, -0.2)],
        "recregions": [fwdpy11.PoissonInterval(0, 1, 1e-3)],
        "gvalue": [fwdpy11.Additive(2.0)],
        "rates": (0, 1e-4, None),
        "simlen": 1000,
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
    }
    return pd


@pytest.mark.parametrize("pop", [{"N": 100, "genome_length": 1}], indirect=["pop"])
def test_single_model_params(pop, pdict1):
    mp = fwdpy11.ModelParams(**pdict1)

    ts = pop.dump_tables_to_tskit(model_params=mp)

    # reconstruct
    mp_rebuilt = fwdpy11.ModelParams(**eval(ts.metadata["model_params"]))

    assert mp == mp_rebuilt


@pytest.mark.parametrize("pop", [{"N": 100, "genome_length": 1}], indirect=["pop"])
def test_multiple_model_params(pop, pdict1, pdict2):
    mp1 = fwdpy11.ModelParams(**pdict1)
    mp2 = fwdpy11.ModelParams(**pdict2)

    ts = pop.dump_tables_to_tskit(model_params={"phase1": mp1, "phase2": mp2})

    # reconstruct
    mp1_rebuilt = fwdpy11.ModelParams(**eval(ts.metadata["model_params"]["phase1"]))
    mp2_rebuilt = fwdpy11.ModelParams(**eval(ts.metadata["model_params"]["phase2"]))

    assert mp1 == mp1_rebuilt
    assert mp2 == mp2_rebuilt


@pytest.mark.parametrize("pop", [{"N": 100, "genome_length": 1}], indirect=["pop"])
def test_demes_graph(pop, gutenkunst):
    ts = pop.dump_tables_to_tskit(demes_graph=gutenkunst)

    # get the graph back out
    g = demes.Graph.fromdict(ts.metadata["demes_graph"])
    assert g == gutenkunst


@pytest.mark.parametrize("pop", [{"N": 100, "genome_length": 1}], indirect=["pop"])
def test_population_metadata(pop, gutenkunst):
    dm = fwdpy11.discrete_demography.from_demes(gutenkunst)
    demes_metadata = {}
    for key, value in dm.metadata["deme_labels"].items():
        demes_metadata[key] = {"name": value}
    ts = pop.dump_tables_to_tskit(population_metadata=demes_metadata)

    for i in range(ts.tables.populations.num_rows):
        assert ts.population(i).metadata == demes_metadata[i]

    # Method 2: get the description out, too
    demes_metadata = {}
    for i, d in enumerate(gutenkunst.demes):
        demes_metadata[i] = {"name": d.name, "description": d.description}

    # Validate that dm thinks that things are in the same order
    for key, value in dm.metadata["deme_labels"].items():
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

    # Test WrappedTreeSequence propery
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning)
        wts = fwdpy11.tskit_tools.WrappedTreeSequence(ts)
    assert wts.data["mydata"] == 11

    class MyType(object):
        def __init__(self, x):
            self.x = x

        def __repr__(self):
            return f"MyType(x={self.x})"

    ts = pop.dump_tables_to_tskit(data=str(MyType(x=11)))
    assert eval(ts.metadata["data"]).x == 11

    # Test WrappedTreeSequence property
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning)
        wts = fwdpy11.tskit_tools.WrappedTreeSequence(ts)
    assert eval(wts.data).x == 11


@pytest.mark.parametrize("pop", [{"N": 100, "genome_length": 1}], indirect=["pop"])
def test_seed(pop):
    ts = pop.dump_tables_to_tskit(seed=333)
    assert ts.metadata["seed"] == 333

    # Test WrappedTreeSequence property
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning)
        wts = fwdpy11.tskit_tools.WrappedTreeSequence(ts)
    assert wts.seed == 333

    with pytest.raises(ValueError):
        _ = pop.dump_tables_to_tskit(seed=-333)
