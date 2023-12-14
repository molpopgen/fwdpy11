import demes
import tskit

from hypothesis import given, settings
from hypothesis.strategies import integers

import fwdpy11


def reconstruct_pop(pop: fwdpy11.DiploidPopulation, ts: tskit.TreeSequence):
    pop2 = fwdpy11.DiploidPopulation.create_from_tskit(
        ts, import_mutations=True, import_individuals=True
    )
    assert pop.generation == pop2.generation
    for i, j in zip(pop.tables.nodes, pop2.tables.nodes):
        assert i.time == j.time
    for i, j in zip(pop.tables.mutations, pop2.tables.mutations):
        assert pop.mutations[i.key] == pop2.mutations[j.key]
        assert pop.mutations[i.key].g == pop2.mutations[j.key].g


@given(seed=integers(0, 42000000), end_time=integers(0, 49))
@settings(deadline=None)
def test_model_that_ends_before_zero(seed, end_time):
    yaml = f"""
    time_units: generations
    demes:
     - name: deme0
       epochs:
        - start_size: 100
          end_time: 50
        - start_size: 100
          end_time: {end_time}
    """
    graph = demes.loads(yaml)
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
    rng = fwdpy11.GSLrng(seed)
    pop = fwdpy11.DiploidPopulation(demog.initial_sizes, 1)
    fwdpy11.evolvets(rng, pop, params, simplification_interval=100)
    ts = pop.dump_tables_to_tskit(model_params=params)
    most_recent_node_time = min(ts.nodes_time)
    assert most_recent_node_time == end_time
    reconstruct_pop(pop, ts)


def test_model_that_ends_before_zero_incompletely_simulated_model():
    yaml = """
    time_units: generations
    demes:
     - name: deme0
       epochs:
        - start_size: 100
          end_time: 50
        - start_size: 100
          end_time: 3
    """
    graph = demes.loads(yaml)
    demog = fwdpy11.ForwardDemesGraph.from_demes(
        graph, burnin=100, burnin_is_exact=True
    )
    pdict = {
        "nregions": [],
        "recregions": [],
        "sregions": [],
        "rates": (0.0, 0.0, 0.0),
        "gvalue": fwdpy11.Multiplicative(scaling=2.0),
        # Only simulate half of the model
        "simlen": demog.final_generation // 2,
        "demography": demog,
    }
    params = fwdpy11.ModelParams(**pdict)
    rng = fwdpy11.GSLrng(101)
    pop = fwdpy11.DiploidPopulation(demog.initial_sizes, 1)
    fwdpy11.evolvets(rng, pop, params, simplification_interval=100)
    ts = pop.dump_tables_to_tskit(model_params=params)
    most_recent_node_time = min(ts.nodes_time)
    assert most_recent_node_time == 3 + (demog.final_generation - pop.generation)
    reconstruct_pop(pop, ts)


@given(seed=integers(0, 42000000), end_time=integers(0, 49))
@settings(deadline=None)
def test_model_that_ends_before_zero_with_mutation(seed, end_time):
    yaml = f"""
    time_units: generations
    demes:
     - name: deme0
       epochs:
        - start_size: 100
          end_time: 50
        - start_size: 100
          end_time: {end_time}
    """
    graph = demes.loads(yaml)
    demog = fwdpy11.ForwardDemesGraph.from_demes(
        graph, burnin=100, burnin_is_exact=True
    )
    pdict = {
        "nregions": [],
        "recregions": [],
        "sregions": [fwdpy11.ExpS(0, 1, 1, mean=-1e-3)],
        "rates": (0.0, 0.01, 0.0),
        "gvalue": fwdpy11.Multiplicative(scaling=2.0),
        "simlen": demog.final_generation,
        "demography": demog,
    }
    params = fwdpy11.ModelParams(**pdict)
    rng = fwdpy11.GSLrng(101)
    pop = fwdpy11.DiploidPopulation(demog.initial_sizes, 1)
    fwdpy11.evolvets(rng, pop, params, simplification_interval=100)
    ts = pop.dump_tables_to_tskit(model_params=params)
    most_recent_node_time = min(ts.nodes_time)
    assert most_recent_node_time == end_time
    reconstruct_pop(pop, ts)
