import demes
import fwdpy11
import pytest


@pytest.fixture
def model_triggering_issue_792():
    yaml = """
description: Two deme model where only one makes it to the end.
time_units: generations
demes:
- name: ancestral
  description: ancestral deme, two epochs
  epochs:
  - {end_time: 50, start_size: 50}
- name: deme1
  description: child 1
  epochs:
  - {start_size: 50, end_size: 50, end_time: 10}
  ancestors: [ancestral]
- name: deme2
  description: child 2
  epochs:
  - {start_size: 50, end_size: 50, end_time: 0}
  ancestors: [ancestral]
"""
    g = demes.loads(yaml)
    return g


def test_issue_792(model_triggering_issue_792):
    """
    Test of GitHub Issue 792
    """
    demog = fwdpy11.discrete_demography.from_demes(model_triggering_issue_792, burnin=1)
    pop = fwdpy11.DiploidPopulation(50, 1)

    pdict = {
        "rates": (0, 0, 0),
        "gvalue": fwdpy11.Multiplicative(2.0),
        "simlen": demog.metadata["total_simulation_length"],
        "demography": demog,
    }

    params = fwdpy11.ModelParams(**pdict)

    rng = fwdpy11.GSLrng(100)

    fwdpy11.evolvets(rng, pop, params, 10)

    _ = pop.dump_tables_to_tskit()
