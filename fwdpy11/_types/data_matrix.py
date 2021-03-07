#
# Copyright (C) 2017-2021 Kevin Thornton <krthornt@uci.edu>
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

from typing import List, Tuple

import numpy as np

from .._fwdpy11 import ll_DataMatrix, StateMatrix


class DataMatrix(ll_DataMatrix):
    """
    Represent a sample from a population in a matrix format.

    There are two possible representations of the data:

    1. As a genotype matrix, where individuals are encoded a 0,1, or 2
    copies of the derived mutation. There is one column per diploid here,
    and one row per variable site.

    2. As a haplotype matrix, with two columns per diploid, and each
    column containing a 0 (ancestral) or 1 (derived) label. Each row
    represents a variable site.

    .. versionchanged:: 0.2.0

        Changed layout to row = variable site.
        Changed to match fwdpp 0.7.0 layout where the neutral
        and selected data are represented as a
        :class:`fwdpy11.StateMatrix`
    """

    def __init__(self, *args):
        super(DataMatrix, self).__init__(*args)

    @property
    def ncol(self) -> int:
        """Sample size of the matrix"""
        return self._ncol

    @property
    def neutral(self) -> StateMatrix:
        """:class:`fwdpy11.StateMatrix` for neutral variants"""
        return self._neutral

    @property
    def neutral_keys(self) -> List[int]:
        """Keys for neutral mutations used to generate matrix."""
        return self._neutral_keys

    @property
    def selected(self) -> StateMatrix:
        """:class:`fwdpy11.StateMatrix` for selected variants"""
        return self._selected

    @property
    def selected_keys(self) -> List[int]:
        """Keys for selected mutations used to generate matrix."""
        return self._selected_keys

    def merge(self) -> Tuple[np.ndarray, List[int]]:
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

        positions = np.concatenate((self.neutral.positions, self.selected.positions))
        sorted_pos_indexes = np.argsort(positions)
        merged = np.concatenate(
            (np.array(self.neutral), np.array(self.selected)), axis=0
        )
        keys = np.concatenate((self.neutral_keys, self.selected_keys))
        return merged[sorted_pos_indexes, :], keys[sorted_pos_indexes]

    def _return_sorted_matrix(self, sm, k):
        pos = sm.positions
        posi = np.argsort(pos)
        return np.array(sm, copy=False)[posi, :], k

    def neutral_matrix(self) -> Tuple[np.ndarray, List[int]]:
        """
        :returns: Neutral mutations and keys
        :rtype: tuple

        The neutral matrix is a :class:`numpy.ndarray`

        .. versionadded:: 0.6.1
        """
        return self._return_sorted_matrix(self.neutral, self.neutral_keys)

    def selected_matrix(self) -> Tuple[np.ndarray, List[int]]:
        """
        :returns: Selected mutations and keys
        :rtype: tuple

        The selected matrix is a :class:`numpy.ndarray`

        .. versionadded:: 0.6.1
        """
        return self._return_sorted_matrix(self.selected, self.selected_keys)
