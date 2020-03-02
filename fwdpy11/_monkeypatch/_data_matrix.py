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

import numpy as np


def _merge(self):
    """
    Merge the neutral and selected data.

    :rtype: tuple
    :returns: a :class:`numpy.ndarray` composed of
              neutral and selected mutations.  Rows (mutations)
              are sorted by position.  The second tuple element
              contains the mutation keys.

    .. versionadded:: 0.6.1
    """
    nneutral = len(self.neutral_keys)
    nselected = len(self.selected_keys)

    if nneutral == 0 and nselected > 0:
        return np.array(self.selected), self.selected_keys
    elif nneutral > 0 and nselected == 0:
        return np.array(self.neutral), self.neutral_keys

    positions = np.concatenate((self.neutral.positions,
                                self.selected.positions))
    sorted_pos_indexes = np.argsort(positions)
    merged = np.concatenate(
        (np.array(self.neutral), np.array(self.selected)), axis=0)
    keys = np.concatenate((self.neutral_keys, self.selected_keys))
    return merged[sorted_pos_indexes, :], keys[sorted_pos_indexes]


def _sorted_neutral(self):
    """
    :returns: Neutral mutations and keys sorted by position
    :rtype: tuple

    The neutral matrix is a :class:`numpy.ndarray`

    .. versionadded:: 0.6.1
    """
    pos = self.neutral.positions
    posi = np.argsort(pos)
    return np.array(self.neutral, copy=False)[posi, :], self.neutral_keys[posi]


def _sorted_selected(self):
    """
    :returns: Selected mutations and keys sorted by position
    :rtype: tuple

    The selected matrix is a :class:`numpy.ndarray`

    .. versionadded:: 0.6.1
    """
    pos = self.selected.positions
    posi = np.argsort(pos)
    return np.array(self.selected, copy=False)[posi, :], \
        self.selected_keys[posi, ]


def _patch_data_matrix(dm):
    dm.merge = _merge
    dm.sorted_neutral = property(_sorted_neutral)
    dm.sorted_selected = property(_sorted_selected)
