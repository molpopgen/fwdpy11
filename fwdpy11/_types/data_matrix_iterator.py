#
# Copyright (C) 2017-2020 Kevin Thornton <krthornt@uci.edu>
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

from typing import List, Tuple, Union

import numpy as np

from .._fwdpy11 import ll_DataMatrixIterator
from .._types import TableCollection


class DataMatrixIterator(ll_DataMatrixIterator):
    """
    Efficient iteration over genomic windows.

    This class allows efficient traversal across multiple
    genomic intervals.  The class is iterable,
    and encapsulates the data fields of
    :class:`fwdpy11.DataMatrix`, which are accessible
    as properties.

    To initialize:

    :param tables: A table collection
    :type tables: :class:`fwdpy11.TableCollection`
    :param samples: A list of samples
    :type samples: list or numpy.ndarray
    :param intervals: The :math:`[start, stop)` positions of each interval
    :type intervals: list[tuple]
    :param neutral: If True, include neutral variants
    :type neutral: bool
    :param selected: If True, include selected variants
    :type selected: bool
    :param fixations: (False) If True, include fixations in the sample
    :type fixations: bool


    .. versionadded:: 0.4.4

    .. versionchanged:: 0.5.0

        Initialization longer requires :class:`fwdpy11.MutationVector`

    """

    def __init__(
        self,
        tables: TableCollection,
        samples: Union[List[int], np.ndarray],
        intervals: List[Tuple[float, float]],
        neutral: bool,
        selected: bool,
        fixations=False,
    ):
        super(DataMatrixIterator, self).__init__(
            tables, samples, intervals, neutral, selected, fixations
        )

    @property
    def neutral(self) -> np.ndarray:
        """Genotypes at neutral variants"""
        return self._neutral

    @property
    def neutral_keys(self) -> np.ndarray:
        """Indexes of the neutral variants."""
        return self._neutral_keys

    @property
    def neutral_positions(self) -> np.ndarray:
        """Positions of neutral variants."""
        return self._neutral_positions

    @property
    def selected(self) -> np.ndarray:
        """Genotypes at selected variants."""
        return self._selected

    @property
    def selected_keys(self) -> np.ndarray:
        """Indexes of selected variants."""
        return self._selected_keys

    @property
    def selected_positions(self) -> np.ndarray:
        """Positions of selected variants."""
        return self._selected_positions

    def __next__(self):
        return self._ll_next()

    def __iter__(self):
        return self._ll_iter()
