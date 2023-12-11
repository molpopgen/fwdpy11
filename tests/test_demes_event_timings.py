from dataclasses import dataclass

import demes
import numpy as np

import fwdpy11


@dataclass
class Parents:
    generation: int
    parents: np.array


@dataclass
class RecordParents:
    parents: list

    def __call__(self, pop, _):
        md = np.array(pop.diploid_metadata, copy=False)
        self.parents.append(Parents(pop.generation, np.copy(md["parents"])))


def run_model(yaml, recorder):
    graph = demes.load(yaml)
    demog = fwdpy11.ForwardDemesGraph.from_demes(
        graph, burnin=100, burnin_is_exact=True
    )
    pdict = {
        "nregions": [],
        "recregions": [],
        "sregions": [],
        "rates": (0.0, 0.0, 0.0),
        "gvalue": fwdpy11.Multiplicative(scaling=2.0),
        "simlen": demog.final_generation,
        "demography": demog,
    }
    params = fwdpy11.ModelParams(**pdict)
    rng = fwdpy11.GSLrng(101)
    pop = fwdpy11.DiploidPopulation(demog.initial_sizes, 1)
    fwdpy11.evolvets(rng, pop, params, simplification_interval=100, recorder=recorder)
    return pop


def test_single_pulse():
    yaml = "tests/demes_event_examples/single_pulse.yaml"
    recorder = RecordParents([])
    pop = run_model(yaml, recorder)
    assert pop.generation == 110
    assert (
        np.all(p.parents.flatten() < 100)
        for p in recorder.parents
        if p.generation == pop.generation - 10
    )
    for p in recorder.parents:
        if p.generation != pop.generation - 10:
            flat = p.parents.flatten()
            deme0 = flat[np.where(flat < 100)]
            deme1 = flat[np.where(flat >= 100)]
            assert len(deme0) == 200, f"{p}"
            assert len(deme1) == 200, f"{p}"


def test_burst_of_migration():
    yaml = "tests/demes_event_examples/burst_of_migration.yaml"
    recorder = RecordParents([])
    pop = run_model(yaml, recorder)
    assert pop.generation == 110
    assert (
        np.all(p.parents.flatten() < 100)
        for p in recorder.parents
        if p.generation > pop.generation - 10 and p.generation <= pop.generation - 5
    )
    for p in recorder.parents:
        if p.generation <= pop.generation - 10 or p.generation > pop.generation - 5:
            flat = p.parents.flatten()
            deme0 = flat[np.where(flat < 100)]
            deme1 = flat[np.where(flat >= 100)]
            assert len(deme0) == 200, f"{p}"
            assert len(deme1) == 200, f"{p}"


def test_deme_existence():
    yaml = "tests/demes_event_examples/deme_existence.yaml"
    recorder = None
    pop = run_model(yaml, recorder)
    assert pop.generation == 150
    ts = pop.dump_tables_to_tskit()
    assert all([i.population == 2 for i in ts.nodes() if i.time < 25.0])
    assert all(
        [i.population == 1 for i in ts.nodes() if i.time >= 25.0 and i.time < 50]
    )
    assert all([i.population == 0 for i in ts.nodes() if i.time >= 50])
