import glob

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
                if all(int(np.rint(i)) > 0 for i in non_integer_deme_sizes):
                    for round_it in [None, True]:
                        _ = fwdpy11.ForwardDemesGraph.from_demes(
                            yaml, burnin=0, round_non_integer_sizes=round_it)
                        _ = fwdpy11.discrete_demography.from_demes(
                            good, burnin=0, round_non_integer_sizes=round_it)
                else:
                    with pytest.raises(ValueError):
                        _ = fwdpy11.ForwardDemesGraph.from_demes(
                            yaml, burnin=0, round_non_integer_sizes=True)
                        _ = fwdpy11.discrete_demography.from_demes(
                            good, burnin=0, round_non_integer_sizes=True)
