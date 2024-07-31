import numpy as np
import pytest

import fwdpy11


def multiplicative_fitness_with_pruning(
    N, diploid, haploid_genomes, mutations, mutation_counts
):
    w = 1.0
    genotypes = {}
    for i in [diploid.first, diploid.second]:
        smutations = memoryview(haploid_genomes[i].smutations)
        for m in smutations:
            if m not in genotypes:
                genotypes[m] = 1
            else:
                genotypes[m] += 1
    for m, c in genotypes.items():
        assert c > 0 and c < 3
        if c == 1:
            e = 1.0 + mutations[m].s * mutations[m].h
        else:
            e = 1.0 + 2.0 * mutations[m].s
        w *= e
    return w


def multiplicative_fitness_without_pruning(
    N, diploid, haploid_genomes, mutations, mutation_counts
):
    w = 1.0
    genotypes = {}
    for i in [diploid.first, diploid.second]:
        smutations = memoryview(haploid_genomes[i].smutations)
        for m in smutations:
            if m not in genotypes:
                genotypes[m] = 1
            else:
                genotypes[m] += 1
    for m, c in genotypes.items():
        assert c > 0 and c < 3
        if c == 1:
            e = 1.0 + mutations[m].s * mutations[m].h
        else:
            e = 1.0 + 2.0 * mutations[m].s
        w *= e
    return w


def test_popgen_model_multiplicative_fitness_with_pruning_and_track_mutations():
    L = 1e6
    N = 1000
    pdict = {
        "sregions": [fwdpy11.ExpS(0, L, 1, mean=0.025)],
        "recregions": [fwdpy11.PoissonInterval(0, L, 0.01)],
        "rates": (0, 1e-3, None),
        "demography": fwdpy11.ForwardDemesGraph.tubes([N], burnin=10),
        "gvalue": fwdpy11.Multiplicative(scaling=2.0),
        "prune_selected": True,
        "simlen": N,
    }
    params = fwdpy11.ModelParams(**pdict)
    pop = fwdpy11.DiploidPopulation(N, L)
    rng = fwdpy11.GSLrng(54321)

    class Recorder:
        def __call__(self, pop, _):
            for i, d in enumerate(pop.diploids):
                w = multiplicative_fitness_with_pruning(
                    pop.N, d, pop.haploid_genomes, pop.mutations, pop.mcounts
                )
                assert np.isclose([w], [pop.diploid_metadata[i].w])

                for g in [d.first, d.second]:
                    for m in pop.haploid_genomes[g].smutations:
                        if pop.mcounts[m] == 2 * pop.N:
                            key = pop.mutations[m].key
                            found = False
                            for f, ft in zip(pop.fixations, pop.fixation_times):
                                if f.key == key:
                                    found = True
                                    assert ft == pop.generation, f"{pop.fixation_times}"
                                    break
                            assert found is True

    class Stop:
        def __call__(self, pop, _):
            return len(pop.fixations) > 1

    fwdpy11.evolvets(
        rng,
        pop,
        params,
        track_mutation_counts=True,
        simplification_interval=100,
        recorder=Recorder(),
        stopping_criterion=Stop(),
    )
    assert len(pop.fixations) > 0


@pytest.mark.parametrize("track_index", [(True, True)])
def test_popgen_model_multiplicative_fitness_simplify_each_generation(
    track_index,
):
    L = 1e6
    N = 1000
    pdict = {
        "sregions": [fwdpy11.ExpS(0, L, 1, mean=0.025)],
        "recregions": [fwdpy11.PoissonInterval(0, L, 0.01)],
        "rates": (0, 1e-3, None),
        "demography": fwdpy11.ForwardDemesGraph.tubes([N], burnin=10),
        "gvalue": fwdpy11.Multiplicative(scaling=2.0),
        "prune_selected": True,
        "simlen": N,
    }
    params = fwdpy11.ModelParams(**pdict)
    pop = fwdpy11.DiploidPopulation(N, L)
    rng = fwdpy11.GSLrng(54321)

    class Recorder:
        def __call__(self, pop, _):
            for i, d in enumerate(pop.diploids):
                # NOTE: not checking fitness calculations
                # here b/c it really slows things down

                for g in [d.first, d.second]:
                    for m in pop.haploid_genomes[g].smutations:
                        if pop.mcounts[m] == 2 * pop.N:
                            key = pop.mutations[m].key
                            found = False
                            for f, ft in zip(pop.fixations, pop.fixation_times):
                                if f.key == key:
                                    found = True
                                    assert ft == pop.generation, f"{pop.fixation_times}"
                                    break
                            assert found is True, f"{track_index}"

    track_mutation_counts, suppress_table_indexing = track_index
    fwdpy11.evolvets(
        rng,
        pop,
        params,
        track_mutation_counts=track_mutation_counts,
        simplification_interval=1,
        recorder=Recorder(),
        stopping_criterion=None,
        suppress_table_indexing=suppress_table_indexing,
    )
    assert len(pop.fixations) > 0


@pytest.mark.parametrize("track_index", [(True, False)])
def test_popgen_model_with_neutral_mutations_multiplicative_fitness_simplify_each_generation(
    track_index,
):
    L = 1e6
    N = 100
    pdict = {
        "sregions": [fwdpy11.ExpS(0, L, 1, mean=0.025)],
        "recregions": [fwdpy11.PoissonInterval(0, L, 0.01)],
        "nregions": [fwdpy11.Region(0, L, 1, label=1)],
        "rates": (1e-2, 1e-3, None),
        "demography": fwdpy11.ForwardDemesGraph.tubes([N], burnin=10),
        "gvalue": fwdpy11.Multiplicative(scaling=2.0),
        "prune_selected": True,
        "simlen": 10 * N,
    }
    params = fwdpy11.ModelParams(**pdict)
    pop = fwdpy11.DiploidPopulation(N, L)
    rng = fwdpy11.GSLrng(54321)

    class Recorder:
        any_neutral = False
        neutral_fixation = False

        def __call__(self, pop, _):
            for i, d in enumerate(pop.diploids):
                # NOTE: not checking fitness calculations
                # here b/c it really slows things down

                for c, m in zip(pop.mcounts, pop.mutations):
                    if m.label == 1:
                        if c == 2 * pop.N:
                            self.neutral_fixation = True
                            key = m.key
                            found = False
                            for f, ft in zip(pop.fixations, pop.fixation_times):
                                if f.key == key:
                                    found = True
                                    assert ft == pop.generation, f"{pop.fixation_times}"
                                    break
                            assert (
                                found is True
                            ), f"failed to find neutral fixation {track_index}"

                for m in pop.tables.mutations:
                    if pop.mutations[m.key].label == 1:
                        self.any_neutral = True
                        assert pop.mcounts[m.key] > 0
                        assert pop.mcounts[m.key] < 2 * N
                        assert pop.mutations[m.key].s == 0.0

                for g in [d.first, d.second]:
                    for m in pop.haploid_genomes[g].smutations:
                        if pop.mcounts[m] == 2 * pop.N:
                            key = pop.mutations[m].key
                            found = False
                            for f, ft in zip(pop.fixations, pop.fixation_times):
                                if f.key == key:
                                    found = True
                                    assert (
                                        ft == pop.generation
                                    ), f"failed to find selected fixation {pop.fixation_times}"
                                    break
                            assert found is True, f"{track_index}"

    track_mutation_counts, suppress_table_indexing = track_index
    recorder = Recorder()
    assert recorder.any_neutral is False
    assert recorder.neutral_fixation is False
    fwdpy11.evolvets(
        rng,
        pop,
        params,
        track_mutation_counts=track_mutation_counts,
        simplification_interval=1,
        recorder=recorder,
        stopping_criterion=None,
        suppress_table_indexing=suppress_table_indexing,
    )
    assert recorder.any_neutral is True
    assert recorder.neutral_fixation is True
    assert len(pop.fixations) > 0
    assert len([i for i in pop.fixations if i.label == 0]) > 0
    assert len([i for i in pop.fixations if i.label == 1]) > 0


def test_popgen_model_multiplicative_fitness_without_pruning_and_track_mutations():
    L = 1e6
    N = 1000
    pdict = {
        "sregions": [fwdpy11.ExpS(0, L, 1, mean=0.025)],
        "recregions": [fwdpy11.PoissonInterval(0, L, 0.01)],
        "rates": (0, 1e-3, None),
        "demography": fwdpy11.ForwardDemesGraph.tubes([N], burnin=10),
        "gvalue": fwdpy11.Multiplicative(scaling=2.0),
        "prune_selected": False,
        "simlen": N,
    }
    params = fwdpy11.ModelParams(**pdict)
    pop = fwdpy11.DiploidPopulation(N, L)
    rng = fwdpy11.GSLrng(54321)

    class Recorder:
        def __call__(self, pop, _):
            for i, d in enumerate(pop.diploids):
                w = multiplicative_fitness_without_pruning(
                    pop.N, d, pop.haploid_genomes, pop.mutations, pop.mcounts
                )
                assert np.isclose([w], [pop.diploid_metadata[i].w])

                for g in [d.first, d.second]:
                    for m in pop.haploid_genomes[g].smutations:
                        if pop.mcounts[m] == 2 * pop.N:
                            key = pop.mutations[m].key
                            found = False
                            for f, ft in zip(pop.fixations, pop.fixation_times):
                                if f.key == key:
                                    found = True
                                    assert ft <= pop.generation, f"{pop.fixation_times}"
                                    break
                            assert found is True

    class Stop:
        def __call__(self, pop, _):
            return len(pop.fixations) > 1

    fwdpy11.evolvets(
        rng,
        pop,
        params,
        track_mutation_counts=True,
        simplification_interval=100,
        recorder=Recorder(),
        stopping_criterion=Stop(),
    )
    assert len(pop.fixations) > 1
    for f in pop.fixations:
        nextant = 0
        found = 0
        for g in pop.haploid_genomes:
            if g.n > 0:
                nextant += 1
                for m in g.smutations:
                    if pop.mutations[m].key == f.key:
                        found += 1
                        break
        assert found == nextant
