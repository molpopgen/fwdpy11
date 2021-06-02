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

from ._flags import *  # NOQA
from .metadata import (DiploidMetadata, decode_individual_metadata,
                       decode_mutation_metadata)
from .trees import WrappedTreeSequence


def load(filename: str):
    """
    Load a tree sequence from a file.

    :param filename: Name of the trees file.
    :type filename: str

    :returns: A tree sequence
    :rtype: :class:`fwdpy11.tskit_tools.WrappedTreeSequence`
    """
    import tskit

    ts = tskit.load(filename)
    return WrappedTreeSequence(ts=ts)
