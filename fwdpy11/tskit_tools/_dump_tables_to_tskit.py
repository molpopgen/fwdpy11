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
import typing

import fwdpy11.tskit_tools
import fwdpy11.tskit_tools.metadata_schema
import numpy as np
import tskit


def _initializePopulationTable(node_view, tc):
    tc.populations.metadata_schema = (
        fwdpy11.tskit_tools.metadata_schema.PopulationMetadata
    )
    for i in sorted(np.unique(node_view["deme"])):
        tc.populations.add_row(metadata={"name": "deme" + str(i)})


def _initializeIndividualTable(self, tc):
    """
    Returns node ID -> individual map
    """
    tc.individuals.metadata_schema = (
        fwdpy11.tskit_tools.metadata_schema.IndividualDiploidMetadata
    )
    # First, alive individuals:
    individal_nodes = {}
    num_ind_nodes = 0
    for i, d in enumerate(self.diploid_metadata):
        individal_nodes[2 * i] = i
        individal_nodes[2 * i + 1] = i
        num_ind_nodes += 1
        tc.individuals.add_row(
            flags=fwdpy11.tskit_tools.INDIVIDUAL_IS_ALIVE,
            metadata=fwdpy11.tskit_tools.metadata_schema.generate_individual_metadata(
                d
            ),
        )

    # Now, preserved nodes
    node_time = np.array(self.tables.nodes, copy=False)["time"]
    for i in self.ancient_sample_metadata:
        assert i.nodes[0] not in individal_nodes, "indivudal record error"
        assert i.nodes[1] not in individal_nodes, "indivudal record error"
        individal_nodes[i.nodes[0]] = num_ind_nodes
        individal_nodes[i.nodes[1]] = num_ind_nodes
        num_ind_nodes += 1
        flag = fwdpy11.tskit_tools.INDIVIDUAL_IS_PRESERVED
        if node_time[i.nodes[0]] == 0.0 and node_time[i.nodes[1]] == 0.0:
            flag |= fwdpy11.tskit_tools.INDIVIDUAL_IS_FIRST_GENERATION
        tc.individuals.add_row(
            flags=flag,
            metadata=fwdpy11.tskit_tools.metadata_schema.generate_individual_metadata(
                i
            ),
        )

    return individal_nodes


def _dump_mutation_site_and_site_tables(self, tc: tskit.TableCollection) -> None:
    mpos = np.array([self.mutations[mr.key].pos for mr in self.tables.mutations])
    ancestral_state = np.zeros(len(self.tables.mutations), dtype=np.int8) + ord("0")
    ancestral_state_offset = np.arange(len(self.tables.mutations) + 1, dtype=np.uint32)
    tc.sites.set_columns(
        position=mpos,
        ancestral_state=ancestral_state,
        ancestral_state_offset=ancestral_state_offset,
    )

    tc.mutations.metadata_schema = (
        fwdpy11.tskit_tools.metadata_schema.determine_mutation_metadata_schema(
            self.mutations
        )
    )
    for m in self.tables.mutations:
        if self.mutations[m.key].g != np.iinfo(np.int32).min:
            origin_time = self.generation - self.mutations[m.key].g
        else:
            origin_time = tc.nodes.time[m.node]

        tc.mutations.add_row(
            site=m.site,
            node=m.node,
            derived_state="1",
            time=origin_time,
            metadata=fwdpy11.tskit_tools.metadata_schema.generate_mutation_metadata(
                m, self.mutations
            ),
        )


def _dump_tables_to_tskit(self, parameters: typing.Optional[typing.Dict] = None):
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

    .. versionchanged:: 0.10.0

        Use tskit metadata schema.
        Mutation time is now stored in the tskit.MutationTable column.
        Origin time of mutations is part of the metadata.
    """
    from fwdpy11 import gsl_version, pybind11_version

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
    for i in np.array(self.ancient_sample_metadata, copy=False)["nodes"].flatten():
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

    _dump_mutation_site_and_site_tables(self, tc)

    tc.provenances.add_row(json.dumps(provenance))
    return tc.tree_sequence()


def _add_dump_tables_to_tskit(cls):
    cls.dump_tables_to_tskit = _dump_tables_to_tskit
