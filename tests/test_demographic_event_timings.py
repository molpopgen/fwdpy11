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

import typing
from dataclasses import dataclass

import fwdpy11
import numpy as np


@dataclass
class DemeSize:
    generation: int
    deme: int
    size: int


@dataclass
class NumMigrants:
    generation: int
    num_migrants: typing.List[int]


@dataclass
class NumSelfers:
    generation: int
    num_selfers: int


class Popsizes(object):
    def __init__(self):
        self.popsizes = []

    def __call__(self, pop, _):
        d = pop.deme_sizes()
        for (i, j) in zip(d[0], d[1]):
            self.popsizes.append(DemeSize(pop.generation, i, j))


class Migrants(object):
    def __init__(self):
        self.migrants = []
        self.popsizes = Popsizes()

    def __call__(self, pop, sampler):
        migrants = [0] * 2
        for i, d in enumerate(pop.diploid_metadata):
            if d.deme == 0:
                if d.parents[0] >= 100:
                    migrants[d.deme] += 1
            else:
                if d.parents[0] < 100:
                    migrants[d.deme] += 1
        self.migrants.append(NumMigrants(pop.generation, migrants))

        self.popsizes(pop, sampler)


class Selfers(object):
    def __init__(self):
        self.selfers = []

    def __call__(self, pop, _):
        selfers = 0
        for d in pop.diploid_metadata:
            if d.parents[0] == d.parents[1]:
                selfers += 1
        self.selfers.append(NumSelfers(pop.generation, selfers))


def setup_sim(popsizes, demography, simlen=None, seed=None, genome_length=None):
    if genome_length is None:
        _genome_len = 1.0
    else:
        _genome_len = genome_length

    if simlen is None:
        _simlen = 5
    else:
        _simlen = simlen

    if seed is None:
        _seed = 100
    else:
        _seed = seed

    pdict = {
        "gvalue": [fwdpy11.Multiplicative(2.0)],
        "rates": (0, 0, 0),
        "simlen": _simlen,
        "demography": demography,
    }

    params = fwdpy11.ModelParams(**pdict)
    pop = fwdpy11.DiploidPopulation(popsizes, _genome_len)
    rng = fwdpy11.GSLrng(_seed)

    return rng, pop, params


# Single events happening in isolation:


def test_single_deme_set_size_time_0():
    demog = fwdpy11.DiscreteDemography(
        set_deme_sizes=[fwdpy11.SetDemeSize(when=0, deme=0, new_size=125)]
    )

    rng, pop, params = setup_sim(100, demog)

    recorder = Popsizes()
    fwdpy11.evolvets(rng, pop, params, 50, recorder=recorder)

    # All pop.generation after 0 (the founders) will have the new size
    # Further, only generations 1- are seen by the sampler
    for i in recorder.popsizes:
        assert i.size == 125


def test_single_deme_set_size_time_simlen():
    demog = fwdpy11.DiscreteDemography(
        set_deme_sizes=[fwdpy11.SetDemeSize(when=5, deme=0, new_size=125)]
    )

    rng, pop, params = setup_sim(100, demog)

    recorder = Popsizes()
    fwdpy11.evolvets(rng, pop, params, 50, recorder=recorder)

    # The size change will not be seen, as it
    # the final generation has not reproduced
    for i in recorder.popsizes:
        assert i.size == 100


def test_single_deme_set_size_time_simlen_then_continue_evolving():
    demog = fwdpy11.DiscreteDemography(
        set_deme_sizes=[fwdpy11.SetDemeSize(when=5, deme=0, new_size=125)]
    )

    rng, pop, params = setup_sim(100, demog)

    recorder = Popsizes()
    fwdpy11.evolvets(rng, pop, params, 50, recorder=recorder)
    assert pop.generation == 5

    # The size change will not be seen, as it
    # the final generation has not reproduced
    for i in recorder.popsizes:
        assert i.size == 100

    # Evolve some more, and the pop size change will
    # "hit" the first daughter generation, which is 6
    fwdpy11.evolvets(rng, pop, params, 50, recorder=recorder)
    assert pop.generation == 10

    for i in recorder.popsizes:
        if i.generation < 6:
            assert i.size == 100
        else:
            assert i.size == 125


def test_single_deme_set_selfing_rate_time_0():
    demog = fwdpy11.DiscreteDemography(
        set_selfing_rates=[fwdpy11.SetSelfingRate(when=0, deme=0, S=1.0)]
    )

    rng, pop, params = setup_sim(100, demog)

    recorder = Selfers()
    fwdpy11.evolvets(rng, pop, params, 50, recorder=recorder)

    for i in recorder.selfers:
        assert i.num_selfers == pop.N


def test_single_deme_set_selfing_rate_time_simlen():
    demog = fwdpy11.DiscreteDemography(
        set_selfing_rates=[fwdpy11.SetSelfingRate(when=5, deme=0, S=1.0)]
    )

    rng, pop, params = setup_sim(100, demog)

    recorder = Selfers()
    fwdpy11.evolvets(rng, pop, params, 50, recorder=recorder)

    assert any([i.num_selfers < pop.N for i in recorder.selfers])


def test_single_deme_set_selfing_rate_time_simlen_then_continue_evolving():
    demog = fwdpy11.DiscreteDemography(
        set_selfing_rates=[fwdpy11.SetSelfingRate(when=5, deme=0, S=1.0)]
    )

    rng, pop, params = setup_sim(100, demog)

    recorder = Selfers()
    fwdpy11.evolvets(rng, pop, params, 50, recorder=recorder)

    assert any([i.num_selfers < pop.N for i in recorder.selfers])

    fwdpy11.evolvets(rng, pop, params, 50, recorder=recorder)

    assert all([i.num_selfers == pop.N for i in recorder.selfers if i.generation > 5])


def test_single_deme_set_growth_rate_time_0():
    demog = fwdpy11.DiscreteDemography(
        set_growth_rates=[fwdpy11.SetExponentialGrowth(when=0, deme=0, G=1.05)]
    )

    rng, pop, params = setup_sim(100, demog)
    ancestral_size = pop.N

    recorder = Popsizes()
    fwdpy11.evolvets(rng, pop, params, 50, recorder=recorder)

    for i in recorder.popsizes:
        assert i.size > ancestral_size


def test_single_deme_set_growth_rate_time_simlen():
    demog = fwdpy11.DiscreteDemography(
        set_growth_rates=[fwdpy11.SetExponentialGrowth(when=5, deme=0, G=1.05)]
    )

    rng, pop, params = setup_sim(100, demog)
    ancestral_size = pop.N

    recorder = Popsizes()
    fwdpy11.evolvets(rng, pop, params, 50, recorder=recorder)

    for i in recorder.popsizes:
        assert i.size == ancestral_size


def test_single_deme_set_growth_rate_time_simlen_then_continue_evolving():
    demog = fwdpy11.DiscreteDemography(
        set_growth_rates=[fwdpy11.SetExponentialGrowth(when=5, deme=0, G=1.05)]
    )

    rng, pop, params = setup_sim(100, demog)
    ancestral_size = pop.N

    recorder = Popsizes()
    fwdpy11.evolvets(rng, pop, params, 50, recorder=recorder)

    for i in recorder.popsizes:
        assert i.size == ancestral_size

    fwdpy11.evolvets(rng, pop, params, 50, recorder=recorder)
    assert all([i.size > ancestral_size for i in recorder.popsizes if i.generation > 5])


def test_two_demes_set_migration_rates_time_0():
    set_migration_rates = [fwdpy11.SetMigrationRates(when=0, deme=0, migrates=[0, 1])]
    set_migration_rates.append(
        fwdpy11.SetMigrationRates(when=0, deme=1, migrates=[1, 0])
    )
    demog = fwdpy11.DiscreteDemography(
        set_migration_rates=set_migration_rates,
        migmatrix=np.identity(2),
    )

    rng, pop, params = setup_sim([100, 100], demog)
    ancestral_sizes = pop.deme_sizes()[1]

    recorder = Migrants()
    fwdpy11.evolvets(rng, pop, params, 50, recorder=recorder)

    for i in recorder.migrants:
        for j, n in enumerate(i.num_migrants):
            assert n == ancestral_sizes[j]


def test_two_demes_set_migration_rates_time_simlen():
    set_migration_rates = [fwdpy11.SetMigrationRates(when=5, deme=0, migrates=[0, 1])]
    set_migration_rates.append(
        fwdpy11.SetMigrationRates(when=5, deme=1, migrates=[1, 0])
    )
    demog = fwdpy11.DiscreteDemography(
        set_migration_rates=set_migration_rates,
        migmatrix=np.identity(2),
    )

    rng, pop, params = setup_sim([100, 100], demog)

    recorder = Migrants()
    fwdpy11.evolvets(rng, pop, params, 50, recorder=recorder)

    for i in recorder.migrants:
        assert all([n == 0 for n in i.num_migrants])


def test_two_demes_set_migration_rates_time_simlen_then_keep_evolving():
    set_migration_rates = [fwdpy11.SetMigrationRates(when=5, deme=0, migrates=[0, 1])]
    set_migration_rates.append(
        fwdpy11.SetMigrationRates(when=5, deme=1, migrates=[1, 0])
    )
    demog = fwdpy11.DiscreteDemography(
        set_migration_rates=set_migration_rates,
        migmatrix=np.identity(2),
    )

    rng, pop, params = setup_sim([100, 100], demog)
    ancestral_sizes = pop.deme_sizes()[1]

    recorder = Migrants()
    fwdpy11.evolvets(rng, pop, params, 50, recorder=recorder)

    for i in recorder.migrants:
        assert all([n == 0 for n in i.num_migrants])

    fwdpy11.evolvets(rng, pop, params, 50, recorder=recorder)
    for i in recorder.migrants:
        if i.generation > 5:
            for j, n in enumerate(i.num_migrants):
                assert n == ancestral_sizes[j]


def test_single_deme_split_by_mass_migration_by_copy_time_0():
    mass_migrations = [
        fwdpy11.copy_individuals(when=0, source=0, destination=1, fraction=0.5)
    ]
    demog = fwdpy11.DiscreteDemography(mass_migrations=mass_migrations)
    rng, pop, params = setup_sim(100, demog)
    recorder = Popsizes()
    fwdpy11.evolvets(rng, pop, params, 50, recorder=recorder)
    for i in recorder.popsizes:
        if i.deme == 0:
            assert i.size == 100
        else:
            assert i.deme == 1
            assert i.size == 50


def test_single_deme_split_by_mass_migration_by_copy_time_simlen():
    mass_migrations = [
        fwdpy11.copy_individuals(when=5, source=0, destination=1, fraction=0.5)
    ]
    demog = fwdpy11.DiscreteDemography(mass_migrations=mass_migrations)
    rng, pop, params = setup_sim(100, demog)
    recorder = Popsizes()
    fwdpy11.evolvets(rng, pop, params, 50, recorder=recorder)
    for i in recorder.popsizes:
        assert i.deme == 0
        assert i.size == 100


def test_single_deme_split_by_mass_migration_by_copy_time_simlen_then_keep_evolving():
    mass_migrations = [
        fwdpy11.copy_individuals(when=5, source=0, destination=1, fraction=0.5)
    ]
    demog = fwdpy11.DiscreteDemography(mass_migrations=mass_migrations)
    rng, pop, params = setup_sim(100, demog)
    recorder = Popsizes()
    fwdpy11.evolvets(rng, pop, params, 50, recorder=recorder)
    for i in recorder.popsizes:
        assert i.deme == 0
        assert i.size == 100

    fwdpy11.evolvets(rng, pop, params, 50, recorder=recorder)
    latest_popsizes = [i for i in recorder.popsizes if i.generation > 5]
    for i in latest_popsizes:
        if i.deme == 0:
            assert i.size == 100
        else:
            assert i.deme == 1
            assert i.size == 50


# Simultaenous events applied to the same deme


def test_single_deme_split_by_mass_migration_by_copy_with_size_change_time_0():
    mass_migrations = [
        fwdpy11.copy_individuals(when=0, source=0, destination=1, fraction=0.5)
    ]
    set_deme_sizes = [fwdpy11.SetDemeSize(when=0, deme=1, new_size=200)]
    demog = fwdpy11.DiscreteDemography(
        mass_migrations=mass_migrations, set_deme_sizes=set_deme_sizes
    )
    rng, pop, params = setup_sim(100, demog)
    recorder = Popsizes()
    fwdpy11.evolvets(rng, pop, params, 50, recorder=recorder)
    for i in recorder.popsizes:
        if i.deme == 0:
            assert i.size == 100
        else:
            assert i.deme == 1
            assert i.size == 200


def test_single_deme_split_by_mass_migration_by_copy_with_size_change_time_simlen():
    mass_migrations = [
        fwdpy11.copy_individuals(when=5, source=0, destination=1, fraction=0.5)
    ]
    set_deme_sizes = [fwdpy11.SetDemeSize(when=5, deme=1, new_size=200)]
    demog = fwdpy11.DiscreteDemography(
        mass_migrations=mass_migrations, set_deme_sizes=set_deme_sizes
    )
    rng, pop, params = setup_sim(100, demog)
    recorder = Popsizes()
    fwdpy11.evolvets(rng, pop, params, 50, recorder=recorder)
    for i in recorder.popsizes:
        assert i.deme == 0
        assert i.size == 100


def test_single_deme_split_by_mass_migration_by_copy_with_size_change_time_simlen_then_keep_evolving():
    mass_migrations = [
        fwdpy11.copy_individuals(when=5, source=0, destination=1, fraction=0.5)
    ]
    set_deme_sizes = [fwdpy11.SetDemeSize(when=5, deme=1, new_size=200)]
    demog = fwdpy11.DiscreteDemography(
        mass_migrations=mass_migrations, set_deme_sizes=set_deme_sizes
    )
    rng, pop, params = setup_sim(100, demog)
    recorder = Popsizes()
    fwdpy11.evolvets(rng, pop, params, 50, recorder=recorder)
    for i in recorder.popsizes:
        assert i.deme == 0
        assert i.size == 100
    fwdpy11.evolvets(rng, pop, params, 50, recorder=recorder)
    latest_popsizes = [i for i in recorder.popsizes if i.generation > 5]
    for i in latest_popsizes:
        if i.deme == 0:
            assert i.size == 100
        else:
            assert i.deme == 1
            assert i.size == 200
