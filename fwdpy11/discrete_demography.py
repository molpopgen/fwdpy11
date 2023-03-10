#
# Copyright (C) 2019 Kevin Thornton <krthornt@uci.edu>
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

import itertools
import typing
import warnings
from collections import defaultdict
from dataclasses import dataclass

import attr
import demes
import demes.demes
import intervaltree
import numpy as np

import fwdpy11

from .class_decorators import (
    attr_add_asblack,
    attr_class_pickle_with_super,
    attr_class_to_from_dict,
    attr_class_to_from_dict_no_recurse,
)

from fwdpy11._types.demographic_model_details import DemographicModelDetails


def from_demes(
    dg: typing.Union[str, demes.Graph], burnin: int = 10, *,
    round_non_integer_sizes: typing.Optional[bool] = None
) -> DemographicModelDetails:
    """
    Build a :class:`fwdpy11.DemographicModelDetails` object using demes.
    The deme graph can either be a demes Graph object or a string as the filepath to a
    demes-specifiend YAML demography.

    :param dg: The demes Graph to convert.
    :type dg: demes.Graph or str
    :param burnin: A factor for how many generations to burn in the simulation.
        For a typical demography with a single root, the number of generations
        of burn in is `burnin` times the root deme's population size. For
        models with multiple root demes joined by migration, that population
        size is determined as the size of the metapopulation.
    :type burnin: int

    :rtype: :class:`fwdpy11.DemographicModelDetails`

    .. versionadded:: 0.14.0
    """
    from ._functions import demography_from_demes

    demog = demography_from_demes(
        dg, burnin, round_non_integer_sizes=round_non_integer_sizes)
    return demog
