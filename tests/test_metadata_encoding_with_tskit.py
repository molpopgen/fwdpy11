#
# Copyright (C) 2017-2020 Kevin Thornton <krthornt@uci.edu>
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

import math
import typing

import attr
import fwdpy11.tskit_tools._dump_tables_to_tskit
import fwdpy11.tskit_tools
import numpy as np
import pytest
import tskit

np.random.seed(54321)

MOCK_NODE_TABLE_DTYPE = [("time", np.float), ("deme", np.int32)]
MOCK_EDGE_TABLE_DTYPE = [
    ("left", np.float),
    ("right", np.float),
    ("parent", np.int32),
    ("child", np.int32),
]


@attr.s(auto_attribs=True, frozen=True, kw_only=True, eq=True)
class MockDiploid(object):
    g: float
    e: float
    w: float
    sex: int
    deme: int
    label: int
    parents: typing.List[int]
    geography: typing.List[float]
    nodes: typing.List[int]


@attr.s(auto_attribs=True, frozen=True, kw_only=True)
class MockMutationRecord(object):
    node: int
    key: int
    site: int
    derived_state: int


@attr.s(auto_attribs=True, frozen=True, kw_only=True)
class MockSite(object):
    position: float
    ancestral_state: int


def empty_mock_node_table() -> np.ndarray:
    return np.array([], dtype=MOCK_NODE_TABLE_DTYPE)


def empty_mock_edge_table() -> np.ndarray:
    return np.array([], dtype=MOCK_EDGE_TABLE_DTYPE)


@attr.s(auto_attribs=True)
class MockTables(object):
    genome_length: float = 1.0
    nodes: np.ndarray = empty_mock_node_table()
    edges: np.ndarray = empty_mock_edge_table()
    mutations: typing.List[MockMutationRecord] = attr.Factory(list)
    sites: typing.List[MockSite] = attr.Factory(list)


@attr.s(auto_attribs=True, kw_only=True)
class MockPopulation(object):
    diploid_metadata: typing.List[MockDiploid] = attr.Factory(list)
    ancient_sample_metadata: typing.List[MockDiploid] = attr.Factory(list)
    mutations: fwdpy11.MutationVector = attr.ib()
    tables: MockTables = attr.ib()

    @mutations.default
    def _init_mutations(self):
        return fwdpy11.MutationVector()

    @tables.default
    def _init_tables(self):
        return MockTables()


def add_diploid_metadata_and_nodes(pop: MockPopulation) -> MockPopulation:
    # Add diploids
    demes = [0, 1, 1, 2]
    node_table_data = []
    for i in range(2):
        node_table_data.append((1, demes[i]))
        node_table_data.append((1, demes[i]))
        pop.diploid_metadata.append(
            MockDiploid(
                g=float(i),
                e=float(i ** 2),
                w=float(i ** 2),
                sex=i,
                label=i,
                deme=demes[i],
                parents=[101, 202],  # nonsense data
                geography=[i, 2 * i, 3 * i],  # nonsense data
                nodes=[2 * i, 2 * i + 1],
            )
        )

    # Ancient samples, also generation 0 == first generation preserved.
    for i in range(2, 4):
        node_table_data.append((0, demes[i]))
        node_table_data.append((0, demes[i]))
        pop.ancient_sample_metadata.append(
            MockDiploid(
                g=float(i),
                e=float(i ** 2),
                w=float(i ** 2),
                sex=i,
                label=i,
                deme=demes[i],
                parents=[101, 202],  # nonsense data
                geography=[i, 2 * i, 3 * i],  # nonsense data
                nodes=[2 * i, 2 * i + 1],
            )
        )
    pop.tables.nodes = np.array(node_table_data, dtype=MOCK_NODE_TABLE_DTYPE)

    return pop


def add_mutation(pop: MockPopulation, node: int, with_vectors: bool):
    pass


@pytest.fixture
def mock_population() -> MockPopulation:
    pop = MockPopulation()

    return pop


@pytest.fixture
def tables() -> tskit.TableCollection:
    tc = tskit.TableCollection(1.0)
    return tc


def test_diploid_metadata(mock_population, tables):
    mock_population = add_diploid_metadata_and_nodes(mock_population)

    fwdpy11.tskit_tools._dump_tables_to_tskit._initializeIndividualTable(
        mock_population, tables
    )

    decoded = fwdpy11.tskit_tools.decode_individual_metadata(tables)

    for i in decoded[: len(mock_population.diploid_metadata)]:
        assert i.alive
        assert not i.preserved
        assert not i.first_generation

    for i in decoded[len(mock_population.diploid_metadata) :]:
        assert not i.alive
        assert i.preserved
        assert i.first_generation

    for i, j in zip(
        decoded,
        mock_population.diploid_metadata + mock_population.ancient_sample_metadata,
    ):
        assert i.g == j.g
        assert i.e == j.e
        assert i.w == j.w
        assert i.sex == j.sex
        assert i.deme == j.deme
        assert i.label == j.label
        assert i.geography == j.geography
        assert i.nodes == j.nodes


def test_mutation_metadata_without_vectors(mock_population, tables):
    mock_population.generation = 1

    mock_population.tables.sites.append(MockSite(position=0.1, ancestral_state=0))

    mock_population.tables.mutations.append(
        MockMutationRecord(node=0, key=0, site=0, derived_state=1)
    )

    mock_population.mutations.append(
        fwdpy11.Mutation(pos=0.1, s=-1e-3, h=1, g=mock_population.generation, label=7)
    )

    fwdpy11.tskit_tools._dump_tables_to_tskit._dump_mutation_site_and_site_tables(
        mock_population, tables
    )
    decoded = fwdpy11.tskit_tools.decode_mutation_metadata(tables)

    assert len(decoded) == len(mock_population.tables.mutations)

    for i, j in zip(decoded, mock_population.mutations):
        assert i.pos == j.pos
        assert i.s == j.s
        assert i.h == j.h
        assert i.g == j.g
        assert i.neutral == j.neutral
        assert i.label == j.label
        assert len(i.esizes) == 0
        assert len(i.heffects) == 0


def test_mutation_metadata_with_vectors(mock_population, tables):
    mock_population.generation = 1

    mock_population.tables.sites.append(MockSite(position=0.1, ancestral_state=0))

    mock_population.tables.mutations.append(
        MockMutationRecord(node=0, key=0, site=0, derived_state=1)
    )

    mock_population.mutations.append(
        fwdpy11.Mutation(
            pos=0.1,
            s=-1e-3,
            h=1,
            g=mock_population.generation,
            label=7,
            esizes=np.arange(4).tolist(),
            heffects=np.array([1.0] * 4).tolist(),
        )
    )

    assert len(mock_population.mutations[0].esizes) == 4

    fwdpy11.tskit_tools._dump_tables_to_tskit._dump_mutation_site_and_site_tables(
        mock_population, tables
    )
    decoded = fwdpy11.tskit_tools.decode_mutation_metadata(tables)

    assert len(decoded) == len(mock_population.tables.mutations)

    for i, j in zip(decoded, mock_population.mutations):
        assert i.pos == j.pos
        assert i.s == j.s
        assert i.h == j.h
        assert i.g == j.g
        assert i.neutral == j.neutral
        assert i.label == j.label
        assert len(i.esizes) == 4
        assert len(i.heffects) == 4
        assert np.array_equal(i.esizes, j.esizes)
        assert np.array_equal(i.heffects, j.heffects)


@pytest.mark.parametrize(
    "msprime_seed", np.random.randint(0, np.iinfo(np.uint32).max, 10).tolist()
)
@pytest.mark.parametrize(
    "fp11_seed", np.random.randint(0, np.iinfo(np.uint32).max, 10).tolist()
)
def test_encoding_neutral_mutations_with_msprime_tree_sequences(
    msprime_seed, fp11_seed
):
    """
    This is a test for GitHub issue 656.
    This test is part of release 0.13.0
    Previous versions of fwdpy11 used uint32_t for Mutation.g.
    The unsigned integer caused overflow issues when adding mutations
    to trees with negative node times.
    Mutation.g was changes to int32_t in 0.13.0 as part of a fix.
    """
    import fwdpy11
    import msprime

    ts = msprime.simulate(500, Ne=1000, random_seed=msprime_seed)
    pop = fwdpy11.DiploidPopulation.create_from_tskit(ts)
    rng = fwdpy11.GSLrng(fp11_seed)
    nmuts = fwdpy11.infinite_sites(rng, pop, 1e-3)
    assert nmuts > 0
    for m in pop.tables.mutations:
        assert pop.mutations[m.key].g <= pop.tables.nodes[m.node].time
    _ = pop.dump_tables_to_tskit()
