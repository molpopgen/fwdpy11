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

import copy

import demes
import fwdpy11
import fwdpy11.conditional_models
import msprime
import numpy as np
import pytest

from fwdpy11_test_utilities import seed_list


def make_pop_with_msprime_ancestry(seed):
    N = 100
    rho = 1000
    L = 1.0
    initial_ts = msprime.sim_ancestry(
        samples=N,
        population_size=2 * N,
        recombination_rate=rho / 4 / N,
        random_seed=seed,
        sequence_length=L,
    )
    pop = fwdpy11.DiploidPopulation.create_from_tskit(initial_ts)
    return pop


@pytest.mark.parametrize("count", [-1, 0, 1])
def test_allele_count(count):
    if count <= 0:
        try:
            _ = fwdpy11.conditional_models.AlleleCount(count)
        except ValueError as _:
            pass
        except Exception as e:
            pytest.fail(f"unexpected exception {e}")
    else:
        _ = fwdpy11.conditional_models.AlleleCount(count)


@pytest.mark.parametrize("r", [(-1, 1), (1, -1), (2, 1)])
def test_bad_allele_count_range(r):
    try:
        _ = fwdpy11.conditional_models.AlleleCountRange(*r)
    except ValueError as _:
        pass
    except Exception as e:
        pytest.fail(f"unexpected exception {e}")


@pytest.mark.parametrize("r", [(0.0, 0.1), (0.1, 1.0)])
def test_bad_frequency_range(r):
    try:
        _ = fwdpy11.conditional_models.FrequencyRange(*r)
    except ValueError as _:
        pass
    except Exception as e:
        pytest.fail(f"unexpected exception {e}")


@pytest.mark.parametrize("msprime_seed", seed_list(151251, 5))
@pytest.mark.parametrize("fp11_seed", seed_list(238145238, 5))
def test_deleterious_mutation_remains_present(msprime_seed, fp11_seed):
    pop = make_pop_with_msprime_ancestry(msprime_seed)
    pdict = {
        "recregions": [fwdpy11.PoissonInterval(0, 1, 5e-2)],
        "gvalue": fwdpy11.Multiplicative(2.0),
        "rates": (0, 0, None),
        "prune_selected": False,
        "simlen": 10,
    }
    params = fwdpy11.ModelParams(**pdict)
    mutation_data = fwdpy11.conditional_models.NewMutationParameters(
        frequency=fwdpy11.conditional_models.AlleleCount(1),
        data=fwdpy11.NewMutationData(effect_size=-1e-2, dominance=1),
        position=fwdpy11.conditional_models.PositionRange(left=0.49, right=0.51),
    )
    rng = fwdpy11.GSLrng(fp11_seed)
    try:
        output = fwdpy11.conditional_models.track_added_mutation(
            rng,
            pop,
            params,
            mutation_data,
            when=0,
            until=10,
        )
        assert output.num_descendant_nodes == 1
        assert output.pop is not None
        assert output.mutation_index is not None
        assert output.mutation_index < len(output.pop.mutations)
        assert output.pop.mutations[output.mutation_index].s == -1e-2
    except fwdpy11.conditional_models.AddMutationFailure as _:
        pass
    except Exception as e:
        pytest.fail(f"unexpected exception: {e}")


@pytest.mark.parametrize("msprime_seed", seed_list(151251, 5))
@pytest.mark.parametrize("fp11_seed", seed_list(238145238, 5))
def test_deleterious_mutation_remains_present_with_final_recording(
    msprime_seed, fp11_seed
):
    pop = make_pop_with_msprime_ancestry(msprime_seed)
    pdict = {
        "recregions": [fwdpy11.PoissonInterval(0, 1, 5e-2)],
        "gvalue": fwdpy11.Multiplicative(2.0),
        "rates": (0, 0, None),
        "prune_selected": False,
        "simlen": 20,
    }
    params = fwdpy11.ModelParams(**pdict)
    mutation_data = fwdpy11.conditional_models.NewMutationParameters(
        frequency=fwdpy11.conditional_models.AlleleCount(1),
        data=fwdpy11.NewMutationData(effect_size=-1e-2, dominance=1),
        position=fwdpy11.conditional_models.PositionRange(left=0.49, right=0.51),
    )
    rng = fwdpy11.GSLrng(fp11_seed)
    try:
        output = fwdpy11.conditional_models.track_added_mutation(
            rng,
            pop,
            params,
            mutation_data,
            when=0,
            until=10,
            sampling_policy=fwdpy11.conditional_models.AncientSamplePolicy.COMPLETION,
        )
        assert output.num_descendant_nodes == 1
        assert output.pop.generation == params.simlen
        assert output.pop is not None
        assert output.mutation_index is not None
        assert output.mutation_index < len(output.pop.mutations)
        assert output.pop.mutations[output.mutation_index].s == -1e-2
        assert len(output.pop.ancient_sample_metadata) == pop.N
        for md in output.pop.ancient_sample_metadata:
            for n in md.nodes:
                assert output.pop.tables.nodes[n].time == 10.0
    except fwdpy11.conditional_models.AddMutationFailure as a:
        pytest.fail(f"failed to add singleton variant: {a}")
    except Exception as e:
        pytest.fail(f"unexpected exception: {e}")


@pytest.mark.parametrize("msprime_seed", seed_list(135123, 5))
@pytest.mark.parametrize("fp11_seed", seed_list(5130125, 5))
@pytest.mark.parametrize("alpha", [1000.0])
def test_sweep_from_new_mutation_using_API(msprime_seed, fp11_seed, alpha):
    pop = make_pop_with_msprime_ancestry(msprime_seed)
    pdict = {
        "recregions": [fwdpy11.PoissonInterval(0, 1, 5e-2)],
        "gvalue": fwdpy11.Multiplicative(2.0),
        "rates": (0, 0, None),
        "prune_selected": False,
        "simlen": 10 * pop.N,
    }
    params = fwdpy11.ModelParams(**pdict)
    mutation_data = fwdpy11.conditional_models.NewMutationParameters(
        frequency=fwdpy11.conditional_models.AlleleCount(1),
        data=fwdpy11.NewMutationData(effect_size=alpha / 2 / pop.N, dominance=1),
        position=fwdpy11.conditional_models.PositionRange(left=0.49, right=0.51),
    )
    rng = fwdpy11.GSLrng(fp11_seed)
    output = fwdpy11.conditional_models.selective_sweep(
        rng,
        pop,
        params,
        mutation_data,
        fwdpy11.conditional_models.GlobalFixation(),
        sampling_policy=fwdpy11.conditional_models.AncientSamplePolicy.COMPLETION,
    )
    assert output.num_descendant_nodes == 1
    assert output.pop.generation > pop.generation
    if output.pop is not None:
        assert output.pop.mcounts[output.mutation_index] == 2 * output.pop.N
        _ = output.pop.dump_tables_to_tskit()
        if output.pop.fixation_times[output.mutation_index] < output.pop.generation:
            for md in output.pop.ancient_sample_metadata:
                for n in md.nodes:
                    assert (
                        output.pop.tables.nodes[n].time
                        == output.pop.fixation_times[output.mutation_index]
                    )
            ti = fwdpy11.TreeIterator(
                output.pop.tables,
                output.pop.ancient_sample_nodes,
                update_samples=True,
            )
            for t in ti:
                for mutation in t.mutations():
                    assert mutation.key == output.mutation_index
                    assert t.leaf_counts(mutation.node) == 2 * pop.N
                    for s in t.samples_below(mutation.node):
                        assert (
                            output.pop.tables.nodes[s].time
                            == output.pop.fixation_times[output.mutation_index]
                        )

        ti = fwdpy11.TreeIterator(output.pop.tables, output.pop.alive_nodes)
        for t in ti:
            for mutation in t.mutations():
                assert mutation.key == output.mutation_index
                assert t.leaf_counts(mutation.node) == 2 * pop.N


@pytest.mark.parametrize("msprime_seed", seed_list(135123, 5))
@pytest.mark.parametrize("fp11_seed", seed_list(5130125, 5))
@pytest.mark.parametrize("alpha", [1000.0])
def test_sweep_from_new_mutation_using_API_exit_when_finished(
    msprime_seed, fp11_seed, alpha
):
    pop = make_pop_with_msprime_ancestry(msprime_seed)
    pdict = {
        "recregions": [fwdpy11.PoissonInterval(0, 1, 5e-2)],
        "gvalue": fwdpy11.Multiplicative(2.0),
        "rates": (0, 0, None),
        "prune_selected": False,
        "simlen": 10 * pop.N,
    }
    params = fwdpy11.ModelParams(**pdict)
    mutation_data = fwdpy11.conditional_models.NewMutationParameters(
        frequency=fwdpy11.conditional_models.AlleleCount(1),
        data=fwdpy11.NewMutationData(effect_size=alpha / 2 / pop.N, dominance=1),
        position=fwdpy11.conditional_models.PositionRange(left=0.49, right=0.51),
    )
    rng = fwdpy11.GSLrng(fp11_seed)
    output = fwdpy11.conditional_models.selective_sweep(
        rng,
        pop,
        params,
        mutation_data,
        fwdpy11.conditional_models.GlobalFixation(),
        sampling_policy=fwdpy11.conditional_models.AncientSamplePolicy.COMPLETION,
        return_when_stopping_condition_met=True,
    )
    assert len(output.pop.fixation_times) == 1
    assert output.pop.generation == output.pop.fixation_times[0]


@pytest.mark.parametrize("msprime_seed", seed_list(135123, 5))
@pytest.mark.parametrize("fp11_seed", seed_list(5130125, 5))
@pytest.mark.parametrize("ndescendants", [2, 7, 10, 23, 12, 100])
@pytest.mark.parametrize("alpha", [1000.0])
def test_sweep_from_standing_variation_using_API(
    msprime_seed, fp11_seed, ndescendants, alpha
):
    pop = make_pop_with_msprime_ancestry(msprime_seed)
    pdict = {
        "recregions": [fwdpy11.PoissonInterval(0, 1, 5e-2)],
        "gvalue": fwdpy11.Multiplicative(2.0),
        "rates": (0, 0, None),
        "prune_selected": False,
        "simlen": 10 * pop.N,
    }
    params = fwdpy11.ModelParams(**pdict)
    mutation_data = fwdpy11.conditional_models.NewMutationParameters(
        frequency=fwdpy11.conditional_models.AlleleCount(ndescendants),
        data=fwdpy11.NewMutationData(effect_size=alpha / 2 / pop.N, dominance=1),
        position=fwdpy11.conditional_models.PositionRange(left=0.49, right=0.51),
    )
    finished = False
    rng = fwdpy11.GSLrng(fp11_seed)
    try:
        output = fwdpy11.conditional_models.selective_sweep(
            rng,
            pop,
            params,
            mutation_data,
            fwdpy11.conditional_models.GlobalFixation(),
        )
        assert output.num_descendant_nodes == ndescendants
        assert output.pop.generation > pop.generation
        if output.pop is not None:
            assert output.pop.mcounts[output.mutation_index] == 2 * output.pop.N
            _ = output.pop.dump_tables_to_tskit()
            finished = True
    except fwdpy11.conditional_models.AddMutationFailure as a:
        pass
    except Exception as e:
        pytest.fail(f"unexpected exception: {e}")


def count_mutation(pop, idx, focal_deme):
    count = 0
    if focal_deme is None:
        count = pop.mcounts[idx]
    else:
        for md in pop.diploid_metadata:
            if md.deme == focal_deme:
                for genome in [
                    pop.diploids[md.label].first,
                    pop.diploids[md.label].second,
                ]:
                    if idx in pop.haploid_genomes[genome].smutations:
                        count += 1

    return count


def base_demes_model():
    yaml = """
    time_units: generations
    demes:
    - name: ancestor
      epochs:
        - {end_time: 100, start_size: 100}"""
    return yaml


def single_deme_model():
    yaml = """
    - name: derived
      ancestors: [ancestor]
      epochs:
        - {end_time: 0, start_size: 100}
    """
    return base_demes_model() + yaml


def two_demes_with_migration():
    yaml = """
    - name: derived0
      ancestors: [ancestor]
      epochs:
        - {end_time: 0, start_size: 50}
    - name: derived1
      ancestors: [ancestor]
      epochs:
        - {end_time: 0, start_size: 50}
    migrations:
      - demes: [derived0, derived1]
        rate: 0.005
    """
    return base_demes_model() + yaml


def two_demes_with_migration_and_growth():
    yaml = """
    - name: derived0
      ancestors: [ancestor]
      epochs:
        - {end_time: 0, start_size: 50, end_size: 60}
    - name: derived1
      ancestors: [ancestor]
      epochs:
        - {end_time: 0, start_size: 50, end_size: 55}
    migrations:
      - demes: [derived0, derived1]
        rate: 0.005
    """
    return base_demes_model() + yaml


def two_demes_with_no_migration():
    yaml = """
    - name: derived0
      ancestors: [ancestor]
      epochs:
        - {end_time: 0, start_size: 50}
    - name: derived1
      ancestors: [ancestor]
      epochs:
        - {end_time: 0, start_size: 50}
    """
    return base_demes_model() + yaml


@pytest.mark.parametrize("fp11_seed", seed_list(5130125, 5))
@pytest.mark.parametrize("demes_yaml", [single_deme_model()])
@pytest.mark.parametrize("alpha", [1000.0])
def test_sweep_from_new_mutation_in_single_deme_using_API(fp11_seed, demes_yaml, alpha):
    g = demes.loads(demes_yaml)
    model = fwdpy11.discrete_demography.from_demes(g, burnin=1)
    isizes = [
        model.metadata["initial_sizes"][d]
        for d in sorted(model.metadata["initial_sizes"].keys())
    ]
    pop = fwdpy11.DiploidPopulation(isizes, 1.0)
    pdict = {
        "recregions": [fwdpy11.PoissonInterval(0, 1, 5e-2)],
        "gvalue": fwdpy11.Multiplicative(2.0),
        "rates": (0, 0, None),
        "demography": model,
        "prune_selected": False,
        "simlen": model.metadata["total_simulation_length"],
    }
    params = fwdpy11.ModelParams(**pdict)
    mutation_data = fwdpy11.conditional_models.NewMutationParameters(
        frequency=fwdpy11.conditional_models.AlleleCount(1),
        data=fwdpy11.NewMutationData(effect_size=alpha / 2 / pop.N, dominance=1),
        position=fwdpy11.conditional_models.PositionRange(left=0.49, right=0.51),
    )
    rng = fwdpy11.GSLrng(fp11_seed)
    output = fwdpy11.conditional_models.selective_sweep(
        rng,
        pop,
        params,
        mutation_data,
        fwdpy11.conditional_models.GlobalFixation(),
        when=model.metadata["burnin_time"] + 50,
    )
    assert output.num_descendant_nodes == 1
    assert output.pop.generation > pop.generation
    if output.pop is not None:
        assert (
            output.pop.mutations[output.mutation_index].g == model.metadata["burnin_time"] + 50
        )
        assert output.pop.mcounts[output.mutation_index] == 2 * output.pop.N
        _ = output.pop.dump_tables_to_tskit()
    else:
        pytest.fail("test failed")


@pytest.mark.parametrize("fp11_seed", seed_list(5130125, 5))
@pytest.mark.parametrize("demes_yaml", [two_demes_with_migration()])
@pytest.mark.parametrize("alpha", [1000.0])
def test_sweep_from_new_mutation_with_demography_using_API(
    fp11_seed, demes_yaml, alpha
):
    g = demes.loads(demes_yaml)
    model = fwdpy11.discrete_demography.from_demes(g, burnin=1)
    isizes = [
        model.metadata["initial_sizes"][d]
        for d in sorted(model.metadata["initial_sizes"].keys())
    ]
    pop = fwdpy11.DiploidPopulation(isizes, 1.0)
    pdict = {
        "recregions": [fwdpy11.PoissonInterval(0, 1, 5e-2)],
        "gvalue": fwdpy11.Multiplicative(2.0),
        "rates": (0, 0, None),
        "demography": model,
        "prune_selected": False,
        "simlen": model.metadata["total_simulation_length"],
    }
    params = fwdpy11.ModelParams(**pdict)
    mutation_data = fwdpy11.conditional_models.NewMutationParameters(
        frequency=fwdpy11.conditional_models.AlleleCount(1),
        data=fwdpy11.NewMutationData(effect_size=alpha / 2 / pop.N, dominance=1),
        position=fwdpy11.conditional_models.PositionRange(left=0.49, right=0.51),
    )
    rng = fwdpy11.GSLrng(fp11_seed)
    output = fwdpy11.conditional_models.selective_sweep(
        rng,
        pop,
        params,
        mutation_data,
        fwdpy11.conditional_models.GlobalFixation(),
        when=model.metadata["burnin_time"] + 50,
    )
    assert output.num_descendant_nodes == 1
    assert output.pop.generation > pop.generation
    assert output.pop.generation == params.simlen
    if output.pop is not None:
        assert (
            output.pop.mutations[output.mutation_index].g == model.metadata["burnin_time"] + 50
        )
        assert output.pop.mcounts[output.mutation_index] == 2 * output.pop.N
        _ = output.pop.dump_tables_to_tskit()
    else:
        pytest.fail("test failed")


@pytest.mark.parametrize("fp11_seed", seed_list(5130125, 1))
@pytest.mark.parametrize("demes_yaml", [two_demes_with_migration()])
@pytest.mark.parametrize("alpha", [1000.0])
def test_origination_deme1_fixation_in_deme_2(fp11_seed, demes_yaml, alpha):
    g = demes.loads(demes_yaml)
    model = fwdpy11.discrete_demography.from_demes(g, burnin=1)
    isizes = [
        model.metadata["initial_sizes"][d]
        for d in sorted(model.metadata["initial_sizes"].keys())
    ]
    pop = fwdpy11.DiploidPopulation(isizes, 1.0)
    pdict = {
        "recregions": [fwdpy11.PoissonInterval(0, 1, 5e-2)],
        "gvalue": fwdpy11.Multiplicative(2.0),
        "demography": model,
        "rates": (0, 0, None),
        "prune_selected": False,
        "simlen": model.metadata["total_simulation_length"],
    }
    params = fwdpy11.ModelParams(**pdict)
    mutation_data = fwdpy11.conditional_models.NewMutationParameters(
        deme=1,
        frequency=fwdpy11.conditional_models.AlleleCount(1),
        data=fwdpy11.NewMutationData(effect_size=alpha / 2 / pop.N, dominance=1),
        position=fwdpy11.conditional_models.PositionRange(left=0.49, right=0.51),
    )
    rng = fwdpy11.GSLrng(fp11_seed)
    output = fwdpy11.conditional_models.selective_sweep(
        rng,
        pop,
        params,
        mutation_data,
        fwdpy11.conditional_models.FocalDemeFixation(deme=2),
        when=model.metadata["burnin_time"] + 50,
    )
    assert output.num_descendant_nodes == 1
    assert output.pop.generation > pop.generation
    if output.pop is not None:
        assert (
            output.pop.mutations[output.mutation_index].g == model.metadata["burnin_time"] + 50
        )
        assert (
            count_mutation(output.pop, output.mutation_index, 2)
            == 2 * output.pop.deme_sizes(as_dict=True)[2]
        )
        _ = output.pop.dump_tables_to_tskit()


@pytest.mark.parametrize("fp11_seed", seed_list(5130125, 1))
@pytest.mark.parametrize("demes_yaml", [two_demes_with_migration_and_growth()])
@pytest.mark.parametrize("alpha", [1000.0])
def test_origination_deme1_fixation_in_deme_2_with_growth(fp11_seed, demes_yaml, alpha):
    g = demes.loads(demes_yaml)
    model = fwdpy11.discrete_demography.from_demes(g, burnin=1)
    isizes = [
        model.metadata["initial_sizes"][d]
        for d in sorted(model.metadata["initial_sizes"].keys())
    ]
    pop = fwdpy11.DiploidPopulation(isizes, 1.0)
    pdict = {
        "recregions": [fwdpy11.PoissonInterval(0, 1, 5e-2)],
        "gvalue": fwdpy11.Multiplicative(2.0),
        "demography": model,
        "rates": (0, 0, None),
        "prune_selected": False,
        "simlen": model.metadata["total_simulation_length"],
    }
    params = fwdpy11.ModelParams(**pdict)
    mutation_data = fwdpy11.conditional_models.NewMutationParameters(
        deme=1,
        frequency=fwdpy11.conditional_models.AlleleCount(1),
        data=fwdpy11.NewMutationData(effect_size=alpha / 2 / pop.N, dominance=1),
        position=fwdpy11.conditional_models.PositionRange(left=0.49, right=0.51),
    )
    rng = fwdpy11.GSLrng(fp11_seed)
    output = fwdpy11.conditional_models.selective_sweep(
        rng,
        pop,
        params,
        mutation_data,
        fwdpy11.conditional_models.FocalDemeFixation(deme=2),
        when=model.metadata["burnin_time"] + 50,
    )
    assert output.num_descendant_nodes == 1
    assert output.pop.generation > pop.generation
    del pop
    if output.pop is not None:
        assert (
            output.pop.mutations[output.mutation_index].g == model.metadata["burnin_time"] + 50
        )

        deme_sizes = output.pop.deme_sizes(as_dict=True)

        for deme, n in zip([1, 2], [60, 55]):
            assert deme_sizes[deme] == n

        assert (
            count_mutation(output.pop, output.mutation_index, 2)
            == 2 * output.pop.deme_sizes(as_dict=True)[2]
        )
        _ = output.pop.dump_tables_to_tskit()
    else:
        pytest.fail("mutation did not fix?")


@pytest.mark.parametrize("fp11_seed", seed_list(5130125, 1))
@pytest.mark.parametrize("demes_yaml", [two_demes_with_no_migration()])
@pytest.mark.parametrize("alpha", [1000.0])
def test_origination_deme2_fixation_in_deme_2_no_migration(
    fp11_seed, demes_yaml, alpha
):
    g = demes.loads(demes_yaml)
    model = fwdpy11.discrete_demography.from_demes(g, burnin=1)
    isizes = [
        model.metadata["initial_sizes"][d]
        for d in sorted(model.metadata["initial_sizes"].keys())
    ]
    pop = fwdpy11.DiploidPopulation(isizes, 1.0)
    pdict = {
        "recregions": [fwdpy11.PoissonInterval(0, 1, 5e-2)],
        "gvalue": fwdpy11.Multiplicative(2.0),
        "demography": model,
        "rates": (0, 0, None),
        "prune_selected": False,
        "simlen": model.metadata["total_simulation_length"],
    }
    params = fwdpy11.ModelParams(**pdict)
    mutation_data = fwdpy11.conditional_models.NewMutationParameters(
        deme=2,
        frequency=fwdpy11.conditional_models.AlleleCount(1),
        data=fwdpy11.NewMutationData(effect_size=alpha / 2 / pop.N, dominance=1),
        position=fwdpy11.conditional_models.PositionRange(left=0.49, right=0.51),
    )
    rng = fwdpy11.GSLrng(fp11_seed)
    output = fwdpy11.conditional_models.selective_sweep(
        rng,
        pop,
        params,
        mutation_data,
        fwdpy11.conditional_models.FocalDemeFixation(deme=2),
        when=model.metadata["burnin_time"] + 50,
    )
    assert output.num_descendant_nodes == 1
    assert output.pop.generation > pop.generation
    if output.pop is not None:
        assert (
            output.pop.mutations[output.mutation_index].g == model.metadata["burnin_time"] + 50
        )
        deme_sizes = output.pop.deme_sizes(as_dict=True)
        assert count_mutation(output.pop, output.mutation_index, 2) == 2 * deme_sizes[2]
        assert count_mutation(output.pop, output.mutation_index, 1) == 0
        _ = output.pop.dump_tables_to_tskit()
