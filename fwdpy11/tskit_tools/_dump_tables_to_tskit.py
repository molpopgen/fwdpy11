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
import tskit  # type: ignore


def _initializePopulationTable(
    node_view, population_metadata: typing.Optional[typing.Dict[int, object]], tc
):
    tc.populations.metadata_schema = (
        fwdpy11.tskit_tools.metadata_schema.PopulationMetadata
    )
    # Prior to 0.15.2, we just iterated over unique values.
    # However, that was incorrect (see GitHub issue 792)
    # because the possibility of demes going extinct could leave
    # nodes in the tree in a way that resulted in population IDs
    # being shifted, giving a LibraryError from tskit.
    for i in range(node_view["deme"].max() + 1):
        if population_metadata is not None and i in population_metadata:
            tc.populations.add_row(metadata=population_metadata[i])
        else:
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
        assert i.nodes[0] not in individal_nodes, "individual record error"
        assert i.nodes[1] not in individal_nodes, "individual record error"
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


def _dump_mutation_site_and_site_tables(
    self, tc: tskit.TableCollection, demographic_model_time_offset: float
) -> None:
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
            origin_time = (
                self.generation
                - self.mutations[m.key].g
                + demographic_model_time_offset
            )
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


def _dump_tables_to_tskit(
    self,
    *,
    model_params=None,
    population_metadata=None,
    data=None,
    seed=None,
    parameters=None,
    destructive=False,
) -> tskit.TreeSequence:
    from .._fwdpy11 import gsl_version, pybind11_version

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
    tc.time_units = "generations"

    tc.metadata_schema = fwdpy11.tskit_tools.metadata_schema.TopLevelMetadata

    # Populate the required fields
    top_level_metadata = {"generation": self.generation}

    if model_params is not None:
        assert isinstance(
            model_params, fwdpy11.ModelParams
        ), "model_params must be instance of fwdpy11.ModelParams"
        top_level_metadata["model_params"] = str(model_params.asdict())

    if data is not None:
        top_level_metadata["data"] = data

    if seed is not None:
        if seed < 0:
            raise ValueError(f"seed must be >=0, got {seed}")
        top_level_metadata["seed"] = seed

    tc.metadata = top_level_metadata

    if destructive is True:
        self._clear_haploid_genomes()

    # We must initialize population and individual
    # tables before we can do anything else.
    # Attempting to set population to anything
    # other than -1 in an tskit.NodeTable will
    # raise an exception if the PopulationTable
    # isn't set up.
    _initializePopulationTable(node_view, population_metadata, tc)
    node_to_individual = _initializeIndividualTable(self, tc)
    individual = [-1 for i in range(len(node_view))]
    for k, v in node_to_individual.items():
        individual[k] = v
    flags = [tskit.NODE_IS_SAMPLE] * 2 * self.N + [0] * (len(node_view) - 2 * self.N)
    # Bug fixed in 0.3.1: add preserved nodes to samples list
    for i in np.array(self.ancient_sample_metadata, copy=False)["nodes"].flatten():
        flags[i] = tskit.NODE_IS_SAMPLE

    demographic_model_time_offset = 0.0
    if model_params is not None:
        try:
            demographic_model_time_offset = float(
                model_params.demography._minimal_end_time_from_demes_graph()
            )
            end_time = model_params.demography.final_generation
        except AttributeError:
            demographic_model_time_offset = float(
                model_params.demography.model._minimal_end_time_from_demes_graph()
            )
            end_time = model_params.demography.model.final_generation
        # Account for demographic models that haven't
        # been simulated all the way through.
        demographic_model_time_offset += end_time - self.generation
    assert demographic_model_time_offset >= 0.0
    tc.nodes.set_columns(
        flags=flags,
        time=node_view["time"] + demographic_model_time_offset,
        population=node_view["deme"],
        individual=individual,
    )
    if destructive is True:
        self._clear_diploid_metadata()
        self._clear_ancient_sample_metadata()
        self.tables._clear_nodes()

    _dump_mutation_site_and_site_tables(self, tc, demographic_model_time_offset)
    if destructive is True:
        self._clear_mutations()
        self.tables._clear_sites()
        self.tables._clear_mutations()

    tc.edges.set_columns(
        left=edge_view["left"],
        right=edge_view["right"],
        parent=edge_view["parent"],
        child=edge_view["child"],
    )
    if destructive is True:
        self.tables._clear_edges()

    tc.provenances.add_row(json.dumps(provenance))

    return tc.tree_sequence()


def _add_dump_tables_to_tskit(cls):
    cls.dump_tables_to_tskit = _dump_tables_to_tskit
