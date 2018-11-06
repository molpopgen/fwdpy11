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
import msprime


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
             'heffects': list(m.heffects)}
        muts.append(str(d).encode('utf-8'))
    return msprime.pack_bytes(muts)


def dump_tables_to_msprime(pop):
    node_view = np.array(pop.tables.nodes, copy=True)
    node_view['time'] -= node_view['time'].max()
    node_view['time'][np.where(node_view['time'] != 0.0)[0]] *= -1.0
    edge_view = np.array(pop.tables.edges, copy=False)
    mut_view = np.array(pop.tables.mutations, copy=False)

    flags = [1]*2*pop.N + [0]*(len(node_view) - 2*pop.N)
    tc = msprime.TableCollection(pop.tables.genome_length())
    tc.nodes.set_columns(flags=flags, time=node_view['time'])
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
