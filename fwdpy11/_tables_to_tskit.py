#
# Copyright (C) 2017 Kevin Thornton <krthornt@uci.edu>
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
import numpy as np
import tskit


# TODO: mutation origin times are recorded forwards in time,
# which is likely confusing, so we omit them
def _generate_mutation_metadata(pop):
    muts = []
    for mr in pop.tables.mutations:
        m = pop.mutations[mr.key]
        d = {'s': m.s,
             'h': m.h,
             # 'g': m.g,
             'label': m.label,
             'esizes': list(m.esizes),
             'heffects': list(m.heffects),
             'neutral': m.neutral}
        muts.append(str(d).encode('utf-8'))
    return tskit.pack_bytes(muts)


def _initializePopulationTable(node_view, tc):
    population_metadata = []
    for i in sorted(np.unique(node_view['population'])):
        md = "deme"+str(i)
        population_metadata.append(md.encode("utf-8"))

    pmd, pmdo = tskit.pack_bytes(population_metadata)
    tc.populations.set_columns(metadata=pmd, metadata_offset=pmdo)


def _generate_individual_metadata(dmd, tc):
    strings = []
    for i in dmd:
        d = {'g': i.g,
             'e': i.e,
             'w': i.w,
             'geography': i.geography,
             'parents': i.parents,
             'sex': i.sex,
             'deme': i.deme,
             'label': i.label}
        strings.append(str(d).encode("utf-8"))
    return strings


def _initializeIndividualTable(pop, tc):
    """
    Returns node ID -> individual map
    """
    # First, alive individuals:
    individal_nodes = {}
    for i in range(pop.N):
        individal_nodes[2*i] = i
        individal_nodes[2*i+1] = i
    metadata_strings = _generate_individual_metadata(pop.diploid_metadata, tc)

    # Now, preserved nodes
    num_ind_nodes = pop.N
    for i in pop.ancient_sample_metadata:
        assert i not in individal_nodes, "indivudal record error"
        individal_nodes[i.nodes[0]] = num_ind_nodes
        individal_nodes[i.nodes[1]] = num_ind_nodes
        num_ind_nodes += 1

    metadata_strings.extend(_generate_individual_metadata(
        pop.ancient_sample_metadata, tc))

    md, mdo = tskit.pack_bytes(metadata_strings)
    flags = [0 for i in range(pop.N+len(pop.ancient_sample_metadata))]
    tc.individuals.set_columns(flags=flags, metadata=md, metadata_offset=mdo)
    return individal_nodes


def dump_tables_to_tskit(pop):
    """
    Converts fwdpy11.TableCollection to an
    tskit.TreeSequence
    """
    node_view = np.array(pop.tables.nodes, copy=True)
    node_view['time'] -= node_view['time'].max()
    node_view['time'][np.where(node_view['time'] != 0.0)[0]] *= -1.0
    edge_view = np.array(pop.tables.edges, copy=False)
    mut_view = np.array(pop.tables.mutations, copy=False)

    tc = tskit.TableCollection(pop.tables.genome_length)

    # We must initialize population and individual
    # tables before we can do anything else.
    # Attempting to set population to anything
    # other than -1 in an tskit.NodeTable will
    # raise an exception if the PopulationTable
    # isn't set up.
    _initializePopulationTable(node_view, tc)
    node_to_individual = _initializeIndividualTable(pop, tc)

    individual = [-1 for i in range(len(node_view))]
    for k, v in node_to_individual.items():
        individual[k] = v
    flags = [1]*2*pop.N + [0]*(len(node_view) - 2*pop.N)
    # Bug fixed in 0.3.1: add preserved nodes to samples list
    for i in pop.tables.preserved_nodes:
        flags[i] = 1
    tc.nodes.set_columns(flags=flags, time=node_view['time'],
                         population=node_view['population'],
                         individual=individual)
    tc.edges.set_columns(left=edge_view['left'],
                         right=edge_view['right'],
                         parent=edge_view['parent'],
                         child=edge_view['child'])

    mpos = np.array([pop.mutations[i].pos for i in mut_view['key']])
    ancestral_state = np.zeros(len(mut_view), dtype=np.int8)+ord('0')
    ancestral_state_offset = np.arange(len(mut_view)+1, dtype=np.uint32)
    tc.sites.set_columns(position=mpos,
                         ancestral_state=ancestral_state,
                         ancestral_state_offset=ancestral_state_offset)

    derived_state = np.zeros(len(mut_view), dtype=np.int8)+ord('1')
    md, mdo = _generate_mutation_metadata(pop)
    tc.mutations.set_columns(site=np.arange(len(mpos), dtype=np.int32),
                             node=mut_view['node'],
                             derived_state=derived_state,
                             derived_state_offset=ancestral_state_offset,
                             metadata=md,
                             metadata_offset=mdo)
    return tc.tree_sequence()
