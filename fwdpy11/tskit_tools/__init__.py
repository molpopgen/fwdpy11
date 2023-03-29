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

"""
This module contains various things to help
interact with ``tskit``.

.. versionadded:: 0.8.0
"""

import typing

import numpy as np
import tskit  # type: ignore

from ._flags import *  # NOQA
from .metadata import (  # NOQA
    DiploidMetadata,
    decode_individual_metadata,
    decode_mutation_metadata,
)


def get_toplevel_metadata(ts: tskit.TreeSequence,
                          name: str) -> typing.Optional[object]:
    """
    Extract top-level metadata from a tree sequence.

    :param ts: A tree sequence
    :type ts: tskit.TreeSequence
    :param name: Metadata field name
    :type name: str
    :returns: The metadata field, if `name` is a valid field, or
              `None` otherwise.
    """
    if name in ts.metadata:
        return ts.metadata[name]
    return None


def iterate_timepoints_with_individuals(
    ts, *, decode_metadata=False
):
    """
    Return an iterator over all unique time points with individuals.

    :param ts: A tree sequence
    :type ts: tskit.TreeSequence
    :param decode_metadata: Flag to decode individual metadata or not
    :type decode_metadata: bool

    For each time point a tuple of (time, nodes, metadata) is yielded.

    :param decode_individual_metadata: Whether to return decoded metadata.
    :type decode_individual_metadata: bool

    If `decode_individual_metadata` is `True`, metadata will be stored in
    a :class:`list` of :class:`fwdpy11.tskit_tools.DiploidMetadata`.
    If `False`, `None` will be yielded.
    """

    nodes_in_individuals = []
    node_times = []
    for individual in ts.individuals():
        for node in individual.nodes:
            nodes_in_individuals.append(node)
            node_times.append(ts.node(node).time)

    node_times = np.array(node_times)
    nodes_in_individuals = np.array(nodes_in_individuals)

    unique_node_times = np.unique(node_times)

    for utime in unique_node_times[::-1]:
        # Get the node tables rows in individuals at this time
        x = np.where(node_times == utime)
        node_table_rows = nodes_in_individuals[x]
        assert np.all(
            np.array([ts.node(i).time for i in node_table_rows]) == utime)

        # Get the individuals
        individuals = np.unique(
            np.array([ts.node(i).individual for i in node_table_rows])
        )
        assert not np.any(individuals == tskit.NULL)

        if decode_metadata is True:
            # now, let's decode the individual metadata for this time slice
            decoded_individual_metadata = decode_individual_metadata(
                ts,
                individuals,
            )
        else:
            decoded_individual_metadata = None
        yield utime, node_table_rows, decoded_individual_metadata
