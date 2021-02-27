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
from typing import Dict, Optional, Iterable

import fwdpy11
import fwdpy11._fwdpy11
import numpy as np


class PopulationMixin(object):
    @property
    def ancient_sample_genetic_values(self) -> np.ndarray:
        """
        Return the genetic values for ancient samples as a readonly 2d
        numpy.ndarray.

        Rows are individuals.  Columns are genetic values.

        ..  versionadded 0.3.0

        """
        return self._ancient_sample_genetic_values

    @property
    def fixations(self) -> Iterable[fwdpy11.Mutation]:
        return self._fixations

    @property
    def fixation_times(self) -> Iterable[int]:
        return self._fixation_times

    @property
    def generation(self) -> int:
        return self._generation

    @property
    def genetic_values(self) -> np.ndarray:
        """
        Return the genetic values as a readonly 2d numpy.ndarray.

        Rows are individuals.  Columns are genetic values.

        ..  versionadded 0.3.0

        """
        return self._genetic_values

    @property
    def haploid_genomes(self) -> Iterable[fwdpy11.HaploidGenome]:
        return self._haploid_genomes

    @property
    def mcounts(self) -> Iterable[int]:
        return self._mcounts

    @property
    def mcounts_ancient_samples(self) -> Iterable[int]:
        return self._mcounts_ancient_samples

    def mutation_indexes(self, pos) -> Optional[np.ndarray]:
        """
        Return indexes in :attr:`fwdpy11.DiploidPopulation.mutations`
        that are mutations at postion `pos`.
        """
        return self._mutation_indexes(pos)

    @property
    def mut_lookup(self) -> Dict:
        return self._mut_lookup

    @property
    def mutations(self) -> Iterable[fwdpy11.Mutation]:
        return self._mutations

    @property
    def mutations_ndarray(self) -> np.ndarray:
        """
        Return readonly numpy.ndarray of mutation data.
        The data are returned as a copy.

        The esizes and heffects fields are not part of the
        array, as their size is not constant and therefore
        they canot be part of a numpy "dtype".

        .. versionadded:: 0.4.0
        """
        return self._mutations_ndarray

    @property
    def N(self):
        return self._N

    def find_fixation_by_key(self, key, offset=0):
        """
        Find a fixation by key.

        :param key: A mutation key. See :func:`fwdpy11.Mutation.key`.
        :type key: tuple
        :param offset: Offset to start search in fixation container.
        :type offset: int

        :rtype: int

        :returns: Index of fixation if found, or -1 otherwise.

        .. versionadded:: 0.2.0

        """
        return self._find_fixation_by_key(key, offset)

    def find_mutation_by_key(self, key, offset):
        """
        Find a mutation by key.

        :param key: A mutation key. See :func:`fwdpy11.Mutation.key`.
        :type key: tuple
        :param offset: Offset to start search in mutation container.
        :type offset: int

        :rtype: int

        :returns: Index of mutation if found, or -1 otherwise.

        .. versionadded:: 0.2.0

        """
        return self._find_mutation_by_key(key, offset)
