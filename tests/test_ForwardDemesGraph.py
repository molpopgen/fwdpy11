import glob
import pickle

import numpy as np
import pytest

import demes
import fwdpy11


def test_invalid_demes_spec_modules():
    for bad in glob.glob("tests/demes-spec/test-cases/invalid/*.yaml"):
        with open(bad, "r") as f:
            yaml = "".join(f.readlines())
            # NOTE: cannot test these YAML with demes.Graph b/c
            # they cannot create a valid Graph.
            with pytest.raises(Exception):
                _ = fwdpy11.ForwardDemesGraph.from_demes(yaml, burnin=0)


def test_valid_demes_spec_modules():
    for good in glob.glob("tests/demes-spec/test-cases/valid/*.yaml"):
        graph = demes.load(good)
        non_integer_deme_sizes = []
        for deme in graph.demes:
            for epoch in deme.epochs:
                for size in [epoch.start_size, epoch.end_size]:
                    if np.isfinite(size) and np.modf(size)[0] != 0.0:
                        non_integer_deme_sizes.append(size)
        if len(non_integer_deme_sizes) == 0:
            with open(good, "r") as f:
                yaml = "".join(f.readlines())
                _ = fwdpy11.ForwardDemesGraph.from_demes(yaml, burnin=0)
                _ = fwdpy11.ForwardDemesGraph.from_demes(
                    demes.loads(yaml), burnin=0)
                with pytest.raises(TypeError):
                    _ = fwdpy11.ForwardDemesGraph.from_demes(
                        yaml, burnin=-1)
        else:
            with open(good, "r") as f:
                yaml = "".join(f.readlines())
                with pytest.raises(ValueError):
                    _ = fwdpy11.ForwardDemesGraph.from_demes(
                        yaml, burnin=0, round_non_integer_sizes=False)
                if all(i > 0.5 for i in non_integer_deme_sizes):
                    for round_it in [True]:
                        _ = fwdpy11.ForwardDemesGraph.from_demes(
                            yaml, burnin=0, round_non_integer_sizes=round_it)
                        _ = fwdpy11.discrete_demography.from_demes(
                            good, burnin=0, round_non_integer_sizes=round_it)
                    for round_it in [None, False]:
                        with pytest.raises(ValueError):
                            _ = fwdpy11.ForwardDemesGraph.from_demes(
                                yaml, burnin=0, round_non_integer_sizes=round_it)
                            _ = fwdpy11.discrete_demography.from_demes(
                                good, burnin=0, round_non_integer_sizes=round_it)
                else:
                    for round_it in [None, False]:
                        with pytest.raises(ValueError):
                            _ = fwdpy11.ForwardDemesGraph.from_demes(
                                yaml, burnin=0, round_non_integer_sizes=round_it)
                            _ = fwdpy11.discrete_demography.from_demes(
                                good, burnin=0, round_non_integer_sizes=round_it)


def test_pickle():
    yaml = """
time_units: generations
demes:
- name: deme1
  start_time: .inf
  epochs:
  - {end_size: 1200.0, end_time: 0, start_size: 1200.0}
- name: deme2
  start_time: .inf
  epochs:
  - {end_size: 100, end_time: 400.0, start_size: 100}
pulses:
- {sources: [deme2], dest: deme1, proportions: [0.09999999999999998], time: 400.0}
migrations: []
    """
    g = fwdpy11.ForwardDemesGraph.from_demes(
        yaml, 93, burnin_is_exact=True, round_non_integer_sizes=True)
    p = pickle.dumps(g)
    up = pickle.loads(p)
    assert g == up
    assert g.graph == up.graph

    # NOTE: we are testing an internal detail here:
    # pickling works via asdict and we don't
    # want to include the demes.Graph in that object.
    d = g.asdict()
    assert "graph" not in d
    assert "yaml" in d
