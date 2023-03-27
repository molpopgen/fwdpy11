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
from dataclasses import dataclass
from typing import Dict, List, Optional, Union

from fwdpy11._types.forward_demes_graph import ForwardDemesGraph
from fwdpy11._types.demographic_model_details import DemographicModelDetails
from .diploid_population import DiploidPopulation


@dataclass
class DemographyDebugger(object):
    """
    Efficiently debug demographic events.

    The class executes a mock simulation
    of demographic events in order to quickly
    detect errors. It also generates a "report"
    of the model details for double-checking.

    If model errors are found, a ``ValueError`` exception
    is raised. The goal is to catch the types of errors that would
    result in an exception during a simulation.

    The class is initialized with the following arguments,
    which may be either positional or keyword:

    :param initial_deme_sizes: The initial sizes of each deme
    :type initial_deme_sizes: list or fwdpy11.DiploidPopulation
    :param events: The demographic model
    :type events: fwdpy11.ForwardDemesGraph or
                  fwdpy11.DemographicModelDetails
    :param simlen: The length of the simulation.  Defaults to `None`.
    :type simlen: int
    :param deme_labels: A map from deme index to a printable name
    :type deme_labels: dict

    .. versionadded:: 0.6.0

    .. versionchanged:: 0.8.1

        Initialization now done with attr.
        A list of initial deme sizes is now accepted.
    """

    initial_deme_sizes: Union[List[int], Dict[int, int], DiploidPopulation]
    events: Union[ForwardDemesGraph, DemographicModelDetails]
    simlen: Optional[int] = None
    deme_labels: Optional[Dict] = None
