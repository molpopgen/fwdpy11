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

from .._fwdpy11 import _data_matrix_from_tables, _make_data_matrix
from .._types import DataMatrix, DiploidPopulation, TableCollection


def make_data_matrix(
    pop: DiploidPopulation,
    samples: Union[List, np.ndarray],
    record_neutral: bool,
    record_selected: bool,
) -> DataMatrix:
    """
    Create a :class:`fwdpy11.DataMatrix` from a table collection.

    :param pop: A population
    :type pop: :class:`fwdpy11.PopulationBase`
    :param samples: A list of sample nodes
    :type samples: list
    :param record_neutral: If True, generate data for neutral variants
    :type record_neutral: bool
    :param record_selected: If True, generate data for selected variants
    :type record_selected: bool

    .. deprecated:: 0.3.0

       Prefer :func:`fwdpy11.data_matrix_from_tables`.

    """
    from warnings import warn

    warn(
        "fwdpy11.make_data_matrix is deprecated."
        " Please use fwdpy11.data_matrix_from_tables.",
        FutureWarning,
    )
    return DataMatrix(_make_data_matrix(pop, samples,
                                        record_neutral, record_selected))


def data_matrix_from_tables(
    tables: TableCollection,
    samples: Union[List, np.ndarray],
    *,
    record_neutral: bool = True,
    record_selected: bool = True,
    include_fixations: bool = False,
    begin: float = 0.0,
    end: Optional[float] = None,
) -> DataMatrix:
    """
    Create a :class:`fwdpy11.DataMatrix` from a table collection.

    :param tables: A TableCollection
    :type tables: fwdpy11.TableCollection
    :param samples: A list of sample nodes
    :type samples: list or :class:`numpy.ndarray`
    :param record_neutral: (True) If True, generate data for neutral variants
    :type record_neutral: bool
    :param record_selected: (True) If True, generate data for selected variants
    :type record_selected: bool
    :param include_selected: (False) Whether to include variants fixed in the sample
    :type include_selected: bool
    :param begin: (0.0) Start of range, inclusive
    :param end: (max float) End of range, exclusive

    :rtype: :class:`fwdpy11.DataMatrix`

    .. versionadded:: 0.3.0

    .. versionchanged:: 0.4.1

           Add begin, end options as floats

    .. versionchanged:: 0.5.0

           No longer requires :class:`fwdpy11.MutationVector` argument

    """

    if end is not None:
        _end = end
    else:
        _end = float(np.finfo(np.float64).max)
    return DataMatrix(
        _data_matrix_from_tables(
            tables,
            samples,
            record_neutral,
            record_selected,
            include_fixations,
            begin,
            _end,
        )
    )
