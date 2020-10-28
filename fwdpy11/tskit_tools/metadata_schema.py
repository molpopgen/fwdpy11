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
import typing

import fwdpy11
import tskit

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
