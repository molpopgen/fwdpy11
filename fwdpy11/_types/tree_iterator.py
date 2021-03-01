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

from typing import Iterable, List, Optional, Union

import numpy as np

from .._fwdpy11 import MutationRecord, Site, ll_TreeIterator
from .table_collection import TableCollection


class TreeIterator(ll_TreeIterator):
    """
    Iterate over the marginal trees in a :class:`fwdpy11.TableCollection`

    Positional arguments for initialization:

    :param tables: A table collection
    :type tables: TableCollection
    :param samples: List of sample nodes
    :type samples: list or numpy.ndarray

    The following are keyword only:

    :type update_samples: If `True`, track list of sample nodes for each tree
    :type update_samples: TableCollection
    :param ancient_samples: List of ancient_samples nodes
    :type ancient_samples: List of sample nodes
    :param begin: Start iteration at this position
    :type begin: float
    :param end: End iteration at this position
    :type end: float or None

    If `end` is `None`, then iteration proceeds over all trees with
    after positions :math:`\\geq \\mathrm{begin}`.

    If `ancient_samples` is not `None`, then
    :meth:`TreeIterator.preserved_leaf_counts` will return the number
    of ancient samples below a node.


    .. versionadded 0.3.0

    .. versionchanged:: 0.4.1

        Add begin, end options as floats for initializing

    """

    def __init__(
        self,
        tables: TableCollection,
        samples: Union[List[int], np.ndarray],
        *,
        update_samples: bool = False,
        ancient_samples: Union[List[int], np.ndarray] = None,
        begin: float = 0.0,
        end: Optional[float] = None
    ):
        _end = end
        if _end is None:
            _end = np.finfo(np.float64).max
        if ancient_samples is None:
            super(TreeIterator, self).__init__(
                tables, samples, update_samples, begin, _end
            )
        else:
            super(TreeIterator, self).__init__(
                tables, samples, ancient_samples, update_samples, begin, _end
            )

    def parent(self, u: int) -> int:
        """
        :param u: A node id
        :type u: int

        :return: Parent id of node `u`
        :rtype: int
        """
        return self._parent(u)

    def leaf_counts(self, u: int) -> int:
        """
        :param u: A node id
        :type u: int

        :return: number of samples (leaves) below node `u`
        :rtype: int
        """
        return self._leaf_counts(u)

    def preserved_leaf_counts(self, u: int) -> int:
        """
        :param u: A node id
        :type u: int

        :return: number of ancient samples below node `u`
        :rtype: int
        """
        return self._preserved_leaf_counts(u)

    def left_sib(self, u: int) -> int:
        """
        :param u: A node id
        :type u: int

        :return: Left sib id of node `u`
        :rtype: int
        """
        return self._left_sib(u)

    def right_sib(self, u: int) -> int:
        """
        :param u: A node id
        :type u: int

        :return: Right sib id of node `u`
        :rtype: int
        """
        return self._right_sib(u)

    def left_child(self, u: int) -> int:
        """
        :param u: A node id
        :type u: int

        :return: Left child id of node `u`
        :rtype: int
        """
        return self._left_child(u)

    def right_child(self, u: int) -> int:
        """
        :param u: A node id
        :type u: int

        :return: Right child id of node `u`
        :rtype: int
        """
        return self._right_child(u)

    def total_time(self) -> float:
        """
        :return: Sum of branch lengths
        :rtype: float
        """
        return self._total_time(self.tables.nodes)

    def nodes(self) -> np.ndarray:
        """
        Return the nodes in the current tree.

        The return order is preorder.

        .. versionadded:: 0.5.1

        """
        return self._nodes()

    def samples(self) -> np.ndarray:
        """
        :return: The samples list
        :rtype: np.ndarray
        """
        return self._samples()

    def samples_below(self, node: int, sort=False) -> np.ndarray:
        """
        Return the list of samples descending from a node.

        :param node: A node id
        :type node: int
        :param sort: (False) Whether or not to sort sample node IDs.
        :type sort: boolean

        .. note::

            Do not store these sample lists without making a "deep"
            copy.  The internal buffer is re-used.
        """
        return self._samples_below(node, sort)

    def sites(self) -> Iterable[Site]:
        """
        Return iterator over all :class:`fwdpy11.Site` objects
        on the current tree.

        .. versionadded:: 0.5.1

        """
        return self._sites()

    def mutations(self) -> Iterable[MutationRecord]:
        """
        Return iterator over all :class:`fwdpy11.MutationRecord` objects
        on the current tree.

        .. versionadded:: 0.5.1

        """
        return self._mutations()

    @property
    def left(self) -> float:
        """Left coordinate of current tree (inclusive)"""
        return self._left

    @property
    def right(self) -> float:
        """Right coordinate of current tree (exclusive)."""
        return self._right

    @property
    def roots(self) -> np.ndarray:
        """Roots of the current tree"""
        return self._roots

    @property
    def sample_size(self) -> int:
        return self._sample_size()

    @property
    def tables(self) -> TableCollection:
        """The underlying TableCollection"""
        return self._tables
