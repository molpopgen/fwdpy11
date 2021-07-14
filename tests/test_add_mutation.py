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

import fwdpy11
import msprime
import numpy as np
import pytest

from fwdpy11_test_utilities import seed_list

# NOTE: this is copied from test/test_tree_sequences.py
# FIXME: this should be a more general fixture?
@pytest.fixture
def set_up_quant_trait_model(simlen=1.0):
    # TODO add neutral variants
    N = 1000
    demography = fwdpy11.DiscreteDemography()
    rho = 2.0
    # theta = 100.
    # nreps = 500
    # mu = theta/(4*N)
    r = rho / (4 * N)
    Opt = fwdpy11.Optimum
    GSSmo = fwdpy11.GSSmo(
        [Opt(when=0, optimum=0.0, VS=1.0), Opt(when=N, optimum=1.0, VS=1.0)]
    )
    a = fwdpy11.Additive(2.0, GSSmo)
    p = {
        "nregions": [],
        "sregions": [fwdpy11.GaussianS(0, 1, 1, 0.25)],
        "recregions": [fwdpy11.Region(0, 1, 1)],
        "rates": (0.0, 0.025, r),
        "gvalue": a,
        "prune_selected": False,
        "demography": demography,
        "simlen": np.rint(simlen * N).astype(int),
    }
    params = fwdpy11.ModelParams(**p)
    rng = fwdpy11.GSLrng(101 * 45 * 110 * 210)
    pop = fwdpy11.DiploidPopulation(N, 1.0)
    return params, rng, pop


def test_new_mutation_data_bad_construction():
    with pytest.raises(ValueError):
        # Mutation is neutral
        fwdpy11.NewMutationData(effect_size=0.0, dominance=1.0)

    with pytest.raises(ValueError):
        fwdpy11.NewMutationData(
            effect_size=np.nan,
            dominance=1.0,
        )

    with pytest.raises(ValueError):
        fwdpy11.NewMutationData(
            effect_size=-1e-3,
            dominance=np.nan,
        )

    with pytest.raises(ValueError):
        fwdpy11.NewMutationData(
            effect_size=0.0,  # 0
            dominance=np.nan,
            esizes=[0.0] * 4,  # all 0
            heffects=[0.5] * 4,
        )

    with pytest.raises(ValueError):
        fwdpy11.NewMutationData(
            effect_size=1e-3,
            dominance=np.nan,
            esizes=[0.0] * 4,  # length 4
            heffects=[0.5] * 5,  # length 5
        )


@pytest.mark.parametrize("ndescendants", [1, 2, 3, 4, 5, 10])
def test_add_mutation_in_single_deme(ndescendants):
    ts = msprime.sim_ancestry(
        50, population_size=50, sequence_length=10, random_seed=612534
    )
    pop = fwdpy11.DiploidPopulation.create_from_tskit(ts)
    data = fwdpy11.NewMutationData(effect_size=-1e-3, dominance=1.0)
    rng = fwdpy11.GSLrng(42)
    key = pop.add_mutation(rng, ndescendants=ndescendants, data=data)
    # Skip over cases where ndescendants could not be satisfied
    if key is not None:
        assert key == 0
        assert pop.mcounts[key] == ndescendants

        count = 0
        for d in pop.diploids:
            for i in [d.first, d.second]:
                assert i >= 0
                assert i < len(pop.haploid_genomes)
                for k in pop.haploid_genomes[i].smutations:
                    if k == key:
                        count += 1
        assert count == pop.mcounts[key]


@pytest.mark.parametrize("ndescendants", [1, 2, 3, 4, 5, 10])
@pytest.mark.parametrize("deme", [0, 1])
def test_add_mutation_in_two_deme_model(deme, ndescendants, small_split_model):
    demog = msprime.Demography.from_demes(small_split_model)
    ts = msprime.sim_ancestry(
        samples={0: 50, 1: 50},
        random_seed=654321,
        sequence_length=1.0,
        demography=demog,
        recombination_rate=1,
        discrete_genome=False,
    )
    pop = fwdpy11.DiploidPopulation.create_from_tskit(ts)
    data = fwdpy11.NewMutationData(effect_size=-1e-3, dominance=1.0)
    rng = fwdpy11.GSLrng(42)
    key = pop.add_mutation(rng, ndescendants=ndescendants, deme=deme, data=data)
    # Skip over cases where ndescendants could not be satisfied
    if key is not None:
        assert key == 0
        assert pop.mcounts[key] == ndescendants

        count = 0
        for di, d in enumerate(pop.diploids):
            for i in [d.first, d.second]:
                assert i >= 0
                assert i < len(pop.haploid_genomes)
                for k in pop.haploid_genomes[i].smutations:
                    if k == key:
                        assert pop.diploid_metadata[di].deme == deme
                        count += 1
        assert count == pop.mcounts[key]


@pytest.mark.parametrize("ndescendants", [3, 4, 5, 10])
@pytest.mark.parametrize("deme", [0, 1])
@pytest.mark.parametrize("window", [(0.3, 0.6), (0.7, 0.88)])
@pytest.mark.parametrize("msprime_seed", [654321, 12345, 111, 222])
@pytest.mark.parametrize("fwdpy11_seed", [123451])
def test_add_mutation_in_two_deme_model_windows(
    deme, ndescendants, small_split_model, window, msprime_seed, fwdpy11_seed
):
    demog = msprime.Demography.from_demes(small_split_model)
    ts = msprime.sim_ancestry(
        samples={0: 50, 1: 50},
        random_seed=msprime_seed,
        sequence_length=1.0,
        demography=demog,
        recombination_rate=1,
        discrete_genome=False,
    )
    pop = fwdpy11.DiploidPopulation.create_from_tskit(ts)
    data = fwdpy11.NewMutationData(effect_size=-1e-3, dominance=1.0)
    rng = fwdpy11.GSLrng(fwdpy11_seed)
    key = pop.add_mutation(
        rng, window=window, ndescendants=ndescendants, deme=deme, data=data
    )
    # Skip over cases where ndescendants could not be satisfied
    if key is not None:
        assert key == 0
        assert pop.mcounts[key] == ndescendants

        pos = pop.mutations[key].pos
        assert pos >= window[0]
        assert pos < window[1]

        count = 0
        for di, d in enumerate(pop.diploids):
            for i in [d.first, d.second]:
                assert i >= 0
                assert i < len(pop.haploid_genomes)
                for k in pop.haploid_genomes[i].smutations:
                    if k == key:
                        assert pop.diploid_metadata[di].deme == deme
                        count += 1
        assert count == pop.mcounts[key]


@pytest.mark.parametrize("ndescendants", [1])
@pytest.mark.parametrize("deme", [0])
def test_add_mutation_in_bad_window(deme, ndescendants, small_split_model):
    demog = msprime.Demography.from_demes(small_split_model)
    ts = msprime.sim_ancestry(
        samples={0: 50, 1: 50},
        random_seed=654321,
        sequence_length=1.0,
        demography=demog,
        recombination_rate=1,
        discrete_genome=False,
    )
    pop = fwdpy11.DiploidPopulation.create_from_tskit(ts)
    data = fwdpy11.NewMutationData(effect_size=0.3, dominance=1.0)
    rng = fwdpy11.GSLrng(1234)
    with pytest.raises(ValueError):
        pop.add_mutation(
            rng, window=(-1, pop.tables.genome_length), ndescendants=1, data=data
        )


def test_adding_to_evolved_pop(set_up_quant_trait_model):
    params, rng, pop = set_up_quant_trait_model
    fwdpy11.evolvets(rng, pop, params, simplification_interval=100)
    assert len(pop.tables.mutations) > 0
    data = fwdpy11.NewMutationData(effect_size=0.3, dominance=1.0)
    key = pop.add_mutation(rng, ndescendants=1, data=data)
    assert key is not None
    assert pop.mcounts[key] == 1
    count = 0
    for d in pop.diploids:
        for i in [d.first, d.second]:
            assert i >= 0
            assert i < len(pop.haploid_genomes)
            for k in pop.haploid_genomes[i].smutations:
                if k == key:
                    count += 1
    assert count == pop.mcounts[key]


def test_attempting_to_add_during_a_simulation(set_up_quant_trait_model):
    params, rng, pop = set_up_quant_trait_model

    class InvalidRecorder(object):
        def __init__(self):
            self.rng = fwdpy11.GSLrng(101)
            self.data = fwdpy11.NewMutationData(effect_size=0.3, dominance=1.0)

        def __call__(self, pop, sampler):
            if pop.generation > 5:
                pop.add_mutation(self.rng, data=self.data)

    with pytest.raises(RuntimeError):
        fwdpy11.evolvets(
            rng, pop, params, simplification_interval=100, recorder=InvalidRecorder()
        )


def test_attempt_to_add_to_pop_without_ancestry():
    rng = fwdpy11.GSLrng(42)
    pop = fwdpy11.DiploidPopulation(100, 1.0)
    data = fwdpy11.NewMutationData(effect_size=0.3, dominance=1.0)
    with pytest.raises(ValueError):
        pop.add_mutation(rng, ndescendants=1, data=data)


@pytest.mark.parametrize("msprime_seed", seed_list(135123, 10))
@pytest.mark.parametrize("fp11_seed", seed_list(5130125, 10))
@pytest.mark.parametrize("ndescendants", [2, 7, 10, 23, 100, 257])
def test_add_mutation_to_random_msprime_output(msprime_seed, fp11_seed, ndescendants):
    initial_ts = msprime.sim_ancestry(
        samples=250,
        population_size=500,
        recombination_rate=1e-3,
        random_seed=msprime_seed,
        sequence_length=1.0,
    )
    pop = fwdpy11.DiploidPopulation.create_from_tskit(initial_ts)
    rng = fwdpy11.GSLrng(fp11_seed)
    data = fwdpy11.NewMutationData(effect_size=1e-3, dominance=1.0)
    idx = pop.add_mutation(
        rng, ndescendants=ndescendants, data=data, window=(0.49, 0.51)
    )
    if idx is not None:
        _ = pop.dump_tables_to_tskit()
