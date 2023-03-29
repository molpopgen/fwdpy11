#
# Copyright (C) 2022 Kevin Thornton <krthornt@uci.edu>
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

import msprime
import numpy as np
import pytest

from hypothesis import given
from hypothesis import settings
from hypothesis.strategies import integers

import fwdpy11


def make_position(sites, edge):
    position = float(np.random.uniform(edge.left, edge.right, 1)[0])
    while position in sites:
        position = float(np.random.uniform(edge.left, edge.right, 1)[0])
    sites.add(position)
    return position, sites


def generate_mutation_origin_time(tables, edge):
    origin_time = int(tables.nodes.time[edge.child] +
                      (tables.nodes.time[edge.parent] -
                       tables.nodes.time[edge.child])/2.0)
    # Sometimes a branch is too short to assign an integer
    # time, so we just bail
    if origin_time >= tables.nodes.time[edge.parent] or \
            origin_time < tables.nodes.time[edge.child]:
        return None

    assert origin_time >= tables.nodes.time[edge.child]
    assert origin_time < tables.nodes.time[edge.parent]
    return origin_time


def make_generic_mutation_metadata(origin_time):
    return {'s': 1e-3, 'h': 1.0, 'origin':
            origin_time,
            'label': 0,
            'neutral': 0,
            # fwdpy11 will have to figure out a real key value later.
            'key': 0}


def make_treeseq_with_one_row_containing_bad_metadata(seed, callback):
    ts = msprime.sim_ancestry(5, population_size=10000, random_seed=seed)

    tables = ts.tables

    tables.mutations.metadata_schema = \
        fwdpy11.tskit_tools.metadata_schema.MutationMetadata
    sites = set()
    for edge in tables.edges:
        position, sites = make_position(sites, edge)
        site = tables.sites.add_row(position, '0')
        origin_time = generate_mutation_origin_time(tables, edge)
        if origin_time is not None:
            md = make_generic_mutation_metadata(origin_time)
            callback(md)
            tables.mutations.add_row(
                site, edge.child, '1', time=origin_time, metadata=md)
            sites.add(position)
            break

    assert tables.mutations.num_rows > 0

    tables.sort()

    ts = tables.tree_sequence()

    return ts

# Failures:
#           recrate=0.00010693073272705097, seed2=1, seed=1,


@given(integers(1, int(2**32 - 1)), integers(1, int(2**32 - 1)))
@settings(deadline=None, max_examples=10)
def test_import_msprime_mutations(seed, seed2):
    ts = msprime.sim_ancestry(5, population_size=10000,
                              recombination_rate=1e-5, sequence_length=1.0,
                              discrete_genome=False, random_seed=seed)
    np.random.seed(seed2)
    tables = ts.tables
    tables.mutations.metadata_schema = \
        fwdpy11.tskit_tools.metadata_schema.MutationMetadata
    sites = set()
    for edge in tables.edges:
        position, sites = make_position(sites, edge)
        origin_time = generate_mutation_origin_time(tables, edge)
        if origin_time is not None:
            md = make_generic_mutation_metadata(origin_time)
            site = tables.sites.add_row(position, '0')
            tables.mutations.add_row(
                site, edge.child, '1', time=origin_time, metadata=md)

    tables.sort()

    ts = tables.tree_sequence()
    pop = fwdpy11.DiploidPopulation.create_from_tskit(
        ts, import_mutations=True)

    assert len(pop.mutations) == ts.num_mutations
    assert len(pop.tables.sites) == ts.num_sites

    for site in ts.sites():
        assert len([i for i in pop.mutations if i.pos == site.position]) == 1
        assert len(
            [i for i in pop.tables.sites if i.position == site.position]) == 1
        assert len(
            [i for i in pop.tables.mutations if
             pop.tables.sites[i.key].position == site.position]) == 1

    # All origin times must be < the time of the mutation time
    for m in pop.tables.mutations:
        node_time = pop.tables.nodes[m.node].time
        assert pop.mutations[m.key].g <= node_time

    # manually tally mutation counts and compare to genomes
    counts = {}
    for g in pop.haploid_genomes:
        for m in g.smutations.tolist() + g.mutations.tolist():
            if m in counts:
                counts[m] += 1
            else:
                counts[m] = 1
    for i, j in counts.items():
        assert pop.mcounts[i] == j

    assert np.array(pop.diploids).flatten().tolist() == [
        (2*i, 2*i+1) for i in range(pop.N)]

    # Finally make sure that which samples have mutations
    # matches up with the tskit input
    for t in ts.trees(sample_lists=True):
        for m in t.mutations():
            samples = [i for i in t.samples(m.node)]
            count = len(samples)
            for i, mut in enumerate(pop.mutations):
                if mut.pos == ts.site(m.site).position:
                    assert count == pop.mcounts[i]

                # now, we brute force get a sample
                # list from our pop w/o traversing the trees
                fp11samples = []
                for dip, md in zip(pop.diploids, pop.diploid_metadata):
                    for key in pop.haploid_genomes[dip.first].smutations:
                        if pop.mutations[key].pos == ts.site(m.site).position:
                            fp11samples.append(md.nodes[0])
                    for key in pop.haploid_genomes[dip.second].smutations:
                        if pop.mutations[key].pos == ts.site(m.site).position:
                            fp11samples.append(md.nodes[1])
                samples = sorted(samples)
                fp11samples = sorted(fp11samples)
                assert samples == fp11samples


@given(integers(1, int(2**32 - 1)), integers(1, int(2**32 - 1)))
@settings(deadline=None, max_examples=10)
def test_import_msprime_mutations_bad_metadata(seed1, seed2):
    np.random.seed(seed1)

    def bad_effect_size(md):
        md['s'] = np.nan

    def bad_dominance(md):
        md['h'] = np.nan

    # The "origin time" in the metadata
    # must equal the mutation's time value

    # NOTE: disabled in 0.19.9 to allow
    # creating of pops from fwdpy11-generated
    # mutations
    # def bad_origin(md):
    #     if md['origin'] > 0.0:
    #         md['origin'] = -1*md['origin']
    #     else:
    #         md['origin'] = -1

    # NOTE: disable in 0.19.8 to
    # fix #1109 but this should be
    # a requirement for mutations added "manually"
    # to tskit.
    # def bad_origin2(md):
    #     if md['origin'] > 0.0:
    #         md['origin'] = -1 #int(0.5*md['origin'])
    #     else:
    #         md['origin'] = -1

    # This will not trigger an error:
    # fwdpy11 sets the neutral field based on if effect_size == 0.0
    # def neutrality_mismatch1(md):
    #     md['neutral'] = 1

    # NOTE: bring back bad_origin and bad_origin2 in 0.20.0 if possible
    for callback in [bad_effect_size, bad_dominance]:
        ts = make_treeseq_with_one_row_containing_bad_metadata(seed2, callback)
        with pytest.raises(ValueError):
            _ = fwdpy11.DiploidPopulation.create_from_tskit(
                ts, import_mutations=True)
