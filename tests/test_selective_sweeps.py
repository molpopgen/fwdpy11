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
import fwdpy11
import msprime
import numpy as np

import pytest
from fwdpy11_test_utilities import seed_list


class MutMonitor(object):
    def __init__(self, index, key):
        self.index = index
        self.key = key

    def __call__(self, pop, _):
        if pop.mcounts[self.index] == 0 or pop.mcounts[self.index] == 2 * pop.N:
            return True
        return False


def run_selective_sweep(
    msprime_seed,
    fp11_seed,
    ndescendants,
    N=500,
    alpha=1000,  # 2Ns
    rho=100,
    L=1.0,
    max_attempts=100,
):
    for initial_ts in msprime.sim_ancestry(
        samples=N,
        population_size=1000,
        recombination_rate=rho / 4 / N,
        random_seed=msprime_seed,
        sequence_length=1.0,
        num_replicates=max_attempts,
    ):
        pop = fwdpy11.DiploidPopulation.create_from_tskit(initial_ts)
        pdict = {
            "recregions": [fwdpy11.PoissonInterval(0, 1, 5e-2)],
            "gvalue": fwdpy11.Multiplicative(2.0),
            "rates": (0, 0, None),
            "prune_selected": False,
            "simlen": 10 * pop.N,
        }
        params = fwdpy11.ModelParams(**pdict)

        rng = fwdpy11.GSLrng(fp11_seed)
        data = fwdpy11.NewMutationData(effect_size=alpha / 2 / N, dominance=1.0)
        idx = pop.add_mutation(
            rng, ndescendants=ndescendants, data=data, window=(0.49, 0.51)
        )
        if idx is None:
            continue

        # Make sure we've chosen a valid time
        assert (
            pop.mutations[idx].g <= pop.tables.nodes[pop.tables.mutations[0].node].time
        ), f"{pop.mutations[idx].g} {pop.tables.nodes[pop.tables.mutations[0].node].time}"
        monitor = MutMonitor(idx, pop.mutations[idx].key)

        fixed = False
        pop_with_fixation = None

        while fixed is not True:
            pcopy = copy.deepcopy(pop)
            assert pcopy.generation == 0
            assert len(pcopy.mutations) == 1

            fwdpy11.evolvets(
                rng,
                pcopy,
                params,
                100,
                track_mutation_counts=True,
                stopping_criterion=monitor,
                suppress_table_indexing=True,
            )

            for f in pcopy.fixations:
                if f.key == monitor.key:
                    pop_with_fixation = pcopy
                    fixed = True

        return pop_with_fixation, idx
    return None, None


@pytest.mark.parametrize("msprime_seed", seed_list(135123, 5))
@pytest.mark.parametrize("fp11_seed", seed_list(5130125, 5))
def test_sweep_from_new_mutation(msprime_seed, fp11_seed):
    pop_with_fixation, idx = run_selective_sweep(msprime_seed, fp11_seed, 1)
    if pop_with_fixation is not None:
        assert pop_with_fixation.mcounts[idx] == 2 * pop_with_fixation.N
        _ = pop_with_fixation.dump_tables_to_tskit()


@pytest.mark.parametrize("msprime_seed", seed_list(135123, 5))
@pytest.mark.parametrize("fp11_seed", seed_list(5130125, 5))
@pytest.mark.parametrize("ndescendants", [2, 7, 10, 23, 100, 257])
def test_sweep_from_standing_variation(msprime_seed, fp11_seed, ndescendants):
    pop_with_fixation, idx = run_selective_sweep(msprime_seed, fp11_seed, ndescendants)
    if pop_with_fixation is not None:
        assert pop_with_fixation.mcounts[idx] == 2 * pop_with_fixation.N
        _ = pop_with_fixation.dump_tables_to_tskit()
