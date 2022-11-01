#
# Copyright (C) 2020 Kevin Thornton <krthornt@uci.edu>
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
import typing

import fwdpy11
import tskit  # type: ignore

TopLevelMetadata = tskit.metadata.MetadataSchema(
    {
        "codec": "json",
        "title": "Top-level metadata for table collection/tree sequence.",
        "type": "object",
        "properties": {
            "generation": {
                "type": "integer",
                "description": "The value of pop.generation at the time data"
                "were exported to tskit",
            },
            "model_params": {
                "type": ["string", "object"],
                "description": "One or more fwdpy11.ModelParams instances.",
            },
            "seed": {
                "type": "integer",
                "description": "Random number seed."
                "This is optional because a random number generator"
                "may be called prior to simulation, thus making the"
                "initial seed not capable of reproducing the simulation",
            },
            "demes_graph": {
                "type": "object",
                "description": "A demographic model specified using demes."
                "This information will be redundant with that stored in model_params,"
                "but it may be useful as it allows reconstruction of the YAML file"
                "from the tree sequence.",
            },
            "data": {
                "type": ["object", "string"],
                "description": "This field is reserved for the user to fill.",
            },
        },
        "required": ["generation"],
    }
)

IndividualDiploidMetadata = tskit.metadata.MetadataSchema(
    {
        "codec": "struct",
        "type": "object",
        "properties": {
            "g": {"type": "number", "binaryFormat": "d"},
            "e": {"type": "number", "binaryFormat": "d"},
            "w": {"type": "number", "binaryFormat": "d"},
            "sex": {"type": "number", "binaryFormat": "i"},
            "deme": {"type": "number", "binaryFormat": "i"},
            "label": {"type": "number", "binaryFormat": "Q"},
            "parents": {
                "type": "array",
                "items": {"type": "number", "binaryFormat": "Q"},
                "arrayLengthFormat": "H",
            },
            "geography": {
                "type": "array",
                "items": {"type": "number", "binaryFormat": "d"},
                "arrayLengthFormat": "H",
            },
            "nodes": {
                "type": "array",
                "items": {"type": "number", "binaryFormat": "i"},
                "arrayLengthFormat": "H",
            },
        },
    }
)

PopulationMetadata = tskit.metadata.MetadataSchema(
    {
        "codec": "json",
        "type": "object",
        "name": "Population metadata",
        "properties": {"name": {"type": "string"}},
    }
)

# NOTE: we could use JSON schema here and
# have the array fields not included in "required".
# However:
# 1. That is much slower.
# 2. This metadata contains floating point values.
#    We don't want to lose precision in case anyone
#    wants to try to recalculate diploid phenotypes/fitnesses

MutationMetadata = tskit.metadata.MetadataSchema(
    {
        "codec": "struct",
        "type": ["object", "null"],
        "name": "Mutation metadata",
        "properties": {
            "s": {"type": "number", "binaryFormat": "d"},
            "h": {"type": "number", "binaryFormat": "d"},
            "origin": {"type": "number", "binaryFormat": "i"},
            "neutral": {"type": "number", "binaryFormat": "?"},
            "label": {"type": "number", "binaryFormat": "H"},
            "key": {"type": "number", "binaryFormat": "Q"},
        },
        "additionalProperties": False,
    }
)

_MutationMetaWithVectorsDict = copy.deepcopy(MutationMetadata.schema)

_MutationMetaWithVectorsDict["name"] = "Mutation metadata with vectors"
_MutationMetaWithVectorsDict["properties"]["esizes"] = {
    "type": "array",
    "items": {"type": "number", "binaryFormat": "d"},
}
_MutationMetaWithVectorsDict["properties"]["heffects"] = {
    "type": "array",
    "items": {"type": "number", "binaryFormat": "d"},
}

MutationMetadataWithVectors = tskit.metadata.MetadataSchema(
    _MutationMetaWithVectorsDict
)


def generate_individual_metadata(
    metadata: fwdpy11._fwdpy11.DiploidMetadata,
) -> typing.Dict:
    d = {
        "g": metadata.g,
        "e": metadata.e,
        "w": metadata.w,
        "geography": [i for i in metadata.geography],
        "parents": [i for i in metadata.parents],
        "nodes": [i for i in metadata.nodes],
        "sex": metadata.sex,
        "deme": metadata.deme,
        "label": metadata.label,
    }
    return d


def determine_mutation_metadata_schema(
    mutations: fwdpy11._fwdpy11.MutationVector,
) -> tskit.metadata.MetadataSchema:
    if len(mutations) == 0 or len(mutations[0].esizes) == 0:
        return MutationMetadata

    return MutationMetadataWithVectors


def generate_mutation_metadata(
    mr: fwdpy11._fwdpy11.MutationRecord, mutations: fwdpy11._fwdpy11.MutationVector
) -> typing.Dict:
    m = mutations[mr.key]
    d = {
        "s": m.s,
        "h": m.h,
        "origin": m.g,
        "label": m.label,
        "neutral": int(m.neutral),
        "key": mr.key,
    }
    if len(m.esizes) > 0:
        d["esizes"] = list(m.esizes)
        d["heffects"] = list(m.heffects)
    return d
