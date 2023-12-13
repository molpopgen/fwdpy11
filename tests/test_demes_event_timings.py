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


# WARNING: all code above is included in the manual.
# Any edits to those lines require that the manual be checked
# to see if it renders correctly.


@dataclass
class DemeSizes:
    generation: int
    sizes: dict


@dataclass
class RecordDemeSizes:
    deme_sizes: list

    def __call__(self, pop: fwdpy11.DiploidPopulation, _):
        self.deme_sizes.append(DemeSizes(pop.generation, pop.deme_sizes(as_dict=True)))


def run_model_return_demography(yaml, recorder):
    pop = run_model(yaml, recorder)
    graph = demes.load(yaml)
    demog = fwdpy11.ForwardDemesGraph.from_demes(
        graph, burnin=100, burnin_is_exact=True
    )
    return pop, demog


def validate_recorder(recorder, demog):
    for i, deme in enumerate(demog.deme_time_intervals()):
        # If a deme is a founder deme, generation 0
        # is not seen by our recorder due to the way
        # we have set up our tests. Therefore, 1 fewer
        # generations is seen.
        if deme.start_time == 0:
            offset = 1
        else:
            offset = 0
        recorded = [d.generation for d in recorder.deme_sizes if i in d.sizes]
        assert len(recorded) > 0
        assert (
            len([j for j in recorded if j >= deme.start_time and j < deme.end_time])
            == deme.end_time - deme.start_time - offset
        )
        for j in range(deme.start_time + offset, deme.end_time):
            assert j in recorded


def validate_deme_existence(yaml):
    recorder = RecordDemeSizes([])
    pop, demog = run_model_return_demography(yaml, recorder)
    validate_recorder(recorder, demog)


def test_single_pulse_deme_sizes():
    yaml = "tests/demes_event_examples/single_pulse.yaml"
    validate_deme_existence(yaml)


def test_burst_of_migration_deme_sizes():
    yaml = "tests/demes_event_examples/burst_of_migration.yaml"
    validate_deme_existence(yaml)


def test_deme_existence_deme_sizes():
    yaml = "tests/demes_event_examples/deme_existence.yaml"
    validate_deme_existence(yaml)


def test_generation_times_deme_sizes():
    yaml = "tests/demes_event_examples/deme_existence_generation_time.yaml"
    validate_deme_existence(yaml)
