#
# Copyright (C) 2021 Kevin Thornton <krthornt@uci.edu>
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

import numpy as np
import tskit

from .._fwdpy11 import Mutation
from .metadata import (DiploidMetadata, decode_individual_metadata,
                       decode_mutation_metadata)


class WrappedTreeSequence(object):
    """
    Encapsulates a :class:`tskit.TreeSequence` as `self.ts`.

    Instance methods give additional functionality that usually
    requires interacting with the metadata.

    Instances of this class may be created using a :class:`tskit.TreeSequence`
    or a call to :func:`fwdpy11.tskit_tools.load`.

    .. versionadded:: 0.15.0
    """

    def __init__(self, ts: tskit.TreeSequence):
        found = False
        for row in ts.provenances():
            record = eval(row.record)
            if "software" in record:
                if record["software"]["name"] == "fwdpy11":
                    found = True
        if not found:
            raise ValueError("this tree sequence was not generated by fwdpy11")

        self.ts = ts

    def timepoints_with_individuals(self, *, decode_metadata=False):
        """
        Return an iterator over all unique time points with individuals.

        For each time point a tuple of (time, nodes, metadata) is yielded.

        :param decode_individual_metadata: Whether to return decoded metadata.
        :type decode_individual_metadata: bool

        If `decode_individual_metadata` is `True`, metadata will be stored in
        a :class:`list` of :class:`fwdpy11.tskit_tools.DiploidMetadata`.
        """
        # Get rows of the node table where the nodes are in individuals
        nodes_in_individuals = np.where(self.ts.tables.nodes.individual != tskit.NULL)[
            0
        ]

        # Get the times
        node_times = self.ts.tables.nodes.time[nodes_in_individuals]

        unique_node_times = np.unique(node_times)

        for utime in unique_node_times[::-1]:
            # Get the node tables rows in individuals at this time
            x = np.where(node_times == utime)
            node_table_rows = nodes_in_individuals[x]
            assert np.all(self.ts.tables.nodes.time[node_table_rows] == utime)

            # Get the individuals
            individuals = np.unique(self.ts.tables.nodes.individual[node_table_rows])
            assert not np.any(individuals == tskit.NULL)

            if decode_metadata is True:
                # now, let's decode the individual metadata for this time slice
                decoded_individual_metadata = decode_individual_metadata(
                    self.ts.tables,
                    individuals,
                )
            else:
                decoded_individual_metadata = None
            yield utime, node_table_rows, decoded_individual_metadata

    def decode_individual_metadata(
        self, rows: typing.Optional[typing.Union[int, slice]] = None
    ) -> typing.List[DiploidMetadata]:
        """
        Decode individual metadata.

        See :func:`fwdpy11.tskit_tools.decode_individual_metadata` for details.
        """
        return decode_individual_metadata(self.ts.tables, rows)

    def decode_mutation_metadata(
        self, rows: typing.Optional[typing.Union[int, slice]] = None
    ) -> typing.List[typing.Optional[Mutation]]:
        """
        Decode mutation metadata.

        See :func:`fwdpy11.tskit_tools.decode_mutation_metadata` for details.
        """
        return decode_mutation_metadata(self.ts.tables, rows)
