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

from typing import List, Optional, Union

import numpy as np

from .._fwdpy11 import MutationRecord, Site, ll_VariantIterator
from .table_collection import TableCollection


class VariantIterator(ll_VariantIterator):
    """
    An iterable class for traversing genotypes in a tree sequence.

    To initialize:

    :param tables: The table collection
    :type tables: fwdpy11.TableCollection
    :param samples: Samples list
    :type samples: list or numpy.ndarray

    The following are `kwarg`-only:

    :param begin: (0.0) First position, inclusive.
    :type begin: float
    :param end: (None) Last position, exclusive.
    :type end: float
    :param include_neutral_variants: (True) Include neutral variants during traversal
    :type include_neutral_variants: bool
    :param include_selected_variants: (True) Include selected variants during traversal
    :type include_selected_variants: bool

    .. versionchanged:: 0.4.1

         Add begin, end options as floats

    .. versionchanged:: 0.4.2

         Add include_neutral_variants and include_selected_variants

    .. versionchanged:: 0.5.0

         No longer requires a :class:`fwdpy11.MutationVector`.

    .. versionchanged:: 0.14.0

         Initialization now requires a TableCollection.
         Some initialization parmeters are now keyword-only.

    """

    def __init__(
        self,
        tables: TableCollection,
        samples: Union[List, np.ndarray],
        *,
        begin: float = 0.0,
        end: Optional[float] = None,
        include_neutral_variants: bool = True,
        include_selected_variants: bool = True
    ):
        if end is not None:
            _end = end
        else:
            _end = np.finfo(np.float64).max

        super(VariantIterator, self).__init__(
            tables,
            samples,
            begin,
            _end,
            include_neutral_variants,
            include_selected_variants,
        )

    @property
    def genotypes(self) -> np.ndarray:
        """Genotype array.  Index order is same as sample input"""
        return self._genotypes

    @property
    def position(self) -> float:
        """Current mutation position"""
        return self._position

    @property
    def records(self) -> List[MutationRecord]:
        """:class:`fwdpy11.MutationRecord` objects for the current Site."""
        return self._records

    @property
    def site(self) -> Site:
        """The current :class:`fwdpy11.Site`"""
        return self._site
