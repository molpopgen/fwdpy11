# Copyright (C) 2023 Kevin Thornton <krthornt@uci.edu>
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
import pytest

import fwdpy11
import numpy as np
import tskit


# NOTE: this is redundant with
# internal checks in the C++ back end.
# However, this code makes use of
# different code paths on the C++ side.
def validate_mutation_table(pop):
    # assert(all([i.time > 0 for i in pop.tables.nodes]))
    ti = fwdpy11.TreeIterator(pop.tables, [i for i in range(2*pop.N)])

    for t in ti:
        for m in t.mutations():
            g = pop.mutations[m.key].g
            nt = pop.tables.nodes[m.node].time
            # Internal mutation  times
            # are meant to be "forwards in time",
            # so a variant must always have an origin
            # time <= that of its node time
            assert g <= nt
            if t.parent(m.node) >= 0:
                pt = pop.tables.nodes[t.parent(m.node)].time
                assert g > pt


def pop_rng_params():
    # Set up parameters
    r = 1e-8
    L = 1e8
    # 1 chromosome, of 1 Morgan length
    assert r * L == 1

    N0 = 1000
    VS = 1.0
    mu = 2.5e-3
    optimum = 0.0

    SD = 0.02
    sregions = fwdpy11.GaussianS(0, L, 1, SD)

    simlen = 3

    pop = fwdpy11.DiploidPopulation(N0, L)

    pdict = {
        "nregions": [],
        "sregions": [sregions],
        "recregions": [fwdpy11.BinomialInterval(0, L, 1)],
        "rates": (0.0, mu, None),
        "gvalue": fwdpy11.Additive(
            scaling=2, gvalue_to_fitness=fwdpy11.GSS(optimum=optimum, VS=VS)
        ),
        "simlen": simlen,
        "prune_selected": False,
    }
    params = fwdpy11.ModelParams(**pdict)

    rng = fwdpy11.GSLrng(424242)

    return pop, rng, params


# Example from Aaron Ragsdale.
# Modified to run much faster
@pytest.fixture
def pop():
    pop, rng, params = pop_rng_params()

    fwdpy11.evolvets(rng, pop, params, 100, suppress_table_indexing=True)

    return pop


# Related to GitHub issue 1109.
def test_recreate_pop_from_own_treeseq(pop):
    gen = pop.generation
    ts = pop.dump_tables_to_tskit()
    nm = ts.num_mutations
    assert nm > 0

    pop = fwdpy11.DiploidPopulation.create_from_tskit(
        ts, import_mutations=True)
    assert pop.generation == gen
    assert len(pop.tables.mutations) == nm
    validate_mutation_table(pop)
    assert(all([pop.mutations[i.key].g >= 0 for i in pop.tables.mutations]))


def test_multiple_export_recreate_evolve_steps():
    pop, rng, params = pop_rng_params()
    for i in range(5):
        fwdpy11.evolvets(rng, pop, params, 100, suppress_table_indexing=True)
        ts = pop.dump_tables_to_tskit()
        pop = fwdpy11.DiploidPopulation.create_from_tskit(
            ts, import_mutations=True)
        assert(all([pop.mutations[i.key].g >= 0 for i in pop.tables.mutations]))
    assert pop.generation == 15


def test_start_pop_from_treeseq_without_mutations(pop):
    ts = pop.dump_tables_to_tskit()
    tables = ts.tables
    nm = len(tables.mutations)
    tables.mutations.clear()
    tables.sites.clear()
    assert len(tables.mutations) == 0
    assert len(tables.sites) == 0
    ts = tables.tree_sequence()

    t = ts.first()
    nodes = [i for i in t.nodes() if t.parent(i) != tskit.NULL and i >
             2000 and (ts.node(t.parent(i)).time - ts.node(i).time) > 1]

    np.random.seed(42)
    rnodes = np.random.choice(nodes, size=nm, replace=False)
    sites = set()
    for n in rnodes:
        time = ts.node(n).time
        ptime = ts.node(t.parent(n)).time
        pos = np.random.uniform(t.interval.left, t.interval.right)
        while pos in sites:
            pos = np.random.uniform(t.interval.left, t.interval.right)
        sites.add(pos)
        site = tables.sites.add_row(pos, ancestral_state='0')
        mtime = np.floor(np.random.uniform(time, ptime))
        assert mtime >= time and mtime < ptime
        md = {'s': 1e-3,
              'h': 1.0,
              'origin': int(mtime),
              'label': np.uint16(0),
              'neutral': 0,
              'key': np.uint64(0),
              }
        tables.mutations.add_row(site, n, '1', time=mtime, metadata=md)
    tables.sort()
    ts = tables.tree_sequence()
    pop = fwdpy11.DiploidPopulation.create_from_tskit(
        ts, import_mutations=True)


def test_complex_start_stop_workflow():
    yaml = """
description: split
time_units: generations
demes:
 - name: A
   epochs:
    - start_size: 10
      end_time: 10
 - name: B
   start_time: 10
   ancestors: [A]
   epochs:
    - start_size: 10
      end_size: 50
 - name: C
   start_time: 10
   ancestors: [A]
   epochs:
    - start_size: 10
      end_size: 50
"""
    demography = fwdpy11.ForwardDemesGraph.from_demes(
        yaml, burnin=10, burnin_is_exact=True)
    pdict = {'recregions': [],
             'nregions': [],
             'sregions': [],
             'rates': (0, 0, 0),
             'gvalue': fwdpy11.Multiplicative(2.),
             'simlen': 11,
             'demography': demography
             }
    params = fwdpy11.ModelParams(**pdict)
    pop = fwdpy11.DiploidPopulation(demography.initial_sizes, 1.0)
    rng = fwdpy11.GSLrng(10)
    fwdpy11.evolvets(rng, pop, params, 10)
    assert len(pop.deme_sizes()) == 2
    assert pop.generation == 11
    ts = pop.dump_tables_to_tskit()
    restored = fwdpy11.DiploidPopulation.create_from_tskit(ts)
    pdict['simlen'] = 100
    params = fwdpy11.ModelParams(**pdict)
    # Evolve rest of model
    fwdpy11.evolvets(rng, restored, params, 10)
    model_time = restored.generation
    assert model_time == 20
    ts = restored.dump_tables_to_tskit()
    restored2 = fwdpy11.DiploidPopulation.create_from_tskit(ts)
    fwdpy11.evolvets(rng, restored2, params, 10)
    assert model_time == restored2.generation
