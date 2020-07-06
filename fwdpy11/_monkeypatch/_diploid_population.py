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

import json

import numpy as np
import tskit

import fwdpy11.tskit_tools


def _alive_nodes(self):
    """
    List of alive nodes corresponding to individuals.

    :rtype: numpy.ndarray

    .. vesionadded:: 0.5.3
    """
    md = np.array(self.diploid_metadata, copy=False)
    return md["nodes"].flatten()


def _get_times(self):
    amd = np.array(self.ancient_sample_metadata, copy=False)
    nodes = np.array(self.tables.nodes, copy=False)
    times = nodes["time"][amd["nodes"][:, 0]]
    utimes = np.unique(times)
    return times, utimes


def _traverse_sample_timepoints(self, include_alive=True):
    """
    Return an iterator over all sample time points.
    The iterator yields time, nodes, and metadata.

    :param include_alive: If True, include currently-alive individuals.

    .. versionadded :: 0.5.1

    .. versionchanged:: 0.5.3

       Monkey-patched into pybind11 class
    """
    amd = np.array(self.ancient_sample_metadata, copy=False)
    amd.flags.writeable = False
    times, utimes = _get_times(self)
    for uti in utimes:
        mdidx = np.where(times == uti)[0]
        sample_nodes = amd["nodes"][mdidx].flatten()
        sample_nodes.flags.writeable = False
        mdslice = amd[mdidx][:]
        mdslice.flags.writeable = False
        yield uti, sample_nodes, mdslice

    if include_alive is True:
        md = np.array(self.diploid_metadata, copy=False)
        nodes = md["nodes"].flatten()
        md.flags.writeable = False
        nodes.flags.writeable = False
        yield self.generation, nodes, md


# NOTE: mutation origin times are recorded forwards in time,
# So we convert them into mutation ages.
def _generate_mutation_metadata(self):
    muts = []
    for mr in self.tables.mutations:
        m = self.mutations[mr.key]
        d = {
            "s": m.s,
            "h": m.h,
            "age": self.generation - m.g + 1,
            "label": m.label,
            "esizes": list(m.esizes),
            "heffects": list(m.heffects),
            "neutral": m.neutral,
            "key": mr.key,
        }
        muts.append(str(d).encode("utf-8"))
    return tskit.pack_bytes(muts)


def _initializePopulationTable(node_view, tc):
    population_metadata = []
    for i in sorted(np.unique(node_view["deme"])):
        md = "deme" + str(i)
        population_metadata.append(md.encode("utf-8"))

    pmd, pmdo = tskit.pack_bytes(population_metadata)
    tc.populations.set_columns(metadata=pmd, metadata_offset=pmdo)


def _generate_individual_metadata(dmd, tc):
    strings = []
    for i in dmd:
        d = {
            "g": i.g,
            "e": i.e,
            "w": i.w,
            "geography": i.geography,
            "parents": i.parents,
            "sex": i.sex,
            "deme": i.deme,
            "label": i.label,
        }
        strings.append(str(d).encode("utf-8"))
    return strings


def _initializeIndividualTable(self, tc):
    """
    Returns node ID -> individual map
    """
    # First, alive individuals:
    individal_nodes = {}
    flags = []
    for i in range(self.N):
        individal_nodes[2 * i] = i
        individal_nodes[2 * i + 1] = i
        flags.append(fwdpy11.tskit_tools.INDIVIDUAL_IS_ALIVE)
    metadata_strings = _generate_individual_metadata(self.diploid_metadata, tc)

    # Now, preserved nodes
    num_ind_nodes = self.N
    for i in self.ancient_sample_metadata:
        assert i not in individal_nodes, "indivudal record error"
        individal_nodes[i.nodes[0]] = num_ind_nodes
        individal_nodes[i.nodes[1]] = num_ind_nodes
        num_ind_nodes += 1
        flag = fwdpy11.tskit_tools.INDIVIDUAL_IS_PRESERVED
        if (
            self.tables.nodes[i.nodes[0]].time == 0.0
            and self.tables.nodes[i.nodes[1]].time == 0.0
        ):
            flag |= fwdpy11.tskit_tools.INDIVIDUAL_IS_FIRST_GENERATION
        flags.append(flag)

    metadata_strings.extend(
        _generate_individual_metadata(self.ancient_sample_metadata, tc)
    )

    md, mdo = tskit.pack_bytes(metadata_strings)
    # flags = [0 for i in range(self.N + len(self.ancient_sample_metadata))]
    tc.individuals.set_columns(flags=flags, metadata=md, metadata_offset=mdo)
    return individal_nodes


def _dump_tables_to_tskit(self, parameters=None):
    """
    Dump the population's TableCollection into
    an tskit TreeSequence

    :param parameters: The simulation parameters for the provenance table.
    :type parameters: None or dict

    :rtype: tskit.TreeSequence

    .. versionchanged:: 0.8.2

        Added `parameters`.
        Generate provenance information for return value.
        The provenance information is validated using
        :func:`tskit.validate_provenance`, which may
        raise an exception.
    """
    from fwdpy11 import pybind11_version, gsl_version

    environment = tskit.provenance.get_environment(
        extra_libs={
            "gsl": {"version": gsl_version()["gsl_version"]},
            "pybind11": {"version": pybind11_version()["pybind11_version"]},
        }
    )

    provenance = {
        "schema_version": "1.0.0",
        "software": {"name": "fwdpy11", "version": f"{fwdpy11.__version__}"},
        "environment": environment,
        "parameters": {},
    }

    if parameters is not None:
        provenance["parameters"] = parameters

    tskit.validate_provenance(provenance)

    node_view = np.array(self.tables.nodes, copy=True)
    node_view["time"] -= node_view["time"].max()
    node_view["time"][np.where(node_view["time"] != 0.0)[0]] *= -1.0
    edge_view = np.array(self.tables.edges, copy=False)
    mut_view = np.array(self.tables.mutations, copy=False)

    tc = tskit.TableCollection(self.tables.genome_length)

    # We must initialize population and individual
    # tables before we can do anything else.
    # Attempting to set population to anything
    # other than -1 in an tskit.NodeTable will
    # raise an exception if the PopulationTable
    # isn't set up.
    _initializePopulationTable(node_view, tc)
    node_to_individual = _initializeIndividualTable(self, tc)

    individual = [-1 for i in range(len(node_view))]
    for k, v in node_to_individual.items():
        individual[k] = v
    flags = [1] * 2 * self.N + [0] * (len(node_view) - 2 * self.N)
    # Bug fixed in 0.3.1: add preserved nodes to samples list
    for i in self.tables.preserved_nodes:
        flags[i] = 1
    tc.nodes.set_columns(
        flags=flags,
        time=node_view["time"],
        population=node_view["deme"],
        individual=individual,
    )
    tc.edges.set_columns(
        left=edge_view["left"],
        right=edge_view["right"],
        parent=edge_view["parent"],
        child=edge_view["child"],
    )

    mpos = np.array([self.mutations[i].pos for i in mut_view["key"]])
    ancestral_state = np.zeros(len(mut_view), dtype=np.int8) + ord("0")
    ancestral_state_offset = np.arange(len(mut_view) + 1, dtype=np.uint32)
    tc.sites.set_columns(
        position=mpos,
        ancestral_state=ancestral_state,
        ancestral_state_offset=ancestral_state_offset,
    )

    derived_state = np.zeros(len(mut_view), dtype=np.int8) + ord("1")
    md, mdo = _generate_mutation_metadata(self)
    tc.mutations.set_columns(
        site=np.arange(len(mpos), dtype=np.int32),
        node=mut_view["node"],
        derived_state=derived_state,
        derived_state_offset=ancestral_state_offset,
        metadata=md,
        metadata_offset=mdo,
    )
    tc.provenances.add_row(json.dumps(provenance))
    return tc.tree_sequence()


def _deme_sizes(self, as_dict=False):
    """
    Return the number of individuals in each deme.

    :param as_dict: If True, return results as a dict.
                    If False, numpy arrays are returned

    The default behavior of this function is equivalent to:

    .. code-block:: python

        md = numpy.array(self.diploid_metadata, copy=False)
        return numpy.unique(md['deme'], return_counts=True)

    .. versionadded:: 0.6.0
    """
    md = np.array(self.diploid_metadata, copy=False)
    deme_sizes = np.unique(md["deme"], return_counts=True)
    if as_dict is False:
        return deme_sizes
    return {i: j for i, j in zip(deme_sizes[0], deme_sizes[1])}


def _patch_diploid_population(d):
    d.alive_nodes = property(_alive_nodes)
    d.sample_timepoints = _traverse_sample_timepoints
    d.dump_tables_to_tskit = _dump_tables_to_tskit
    d.deme_sizes = _deme_sizes
