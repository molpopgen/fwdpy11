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


@attr.s(repr_ns="fwdpy11")
@attr_class_pickle_with_super
@attr_class_to_from_dict
class ForwardDemesGraph(fwdpy11._fwdpy11._ForwardDemesGraph):
    """
    A forward-in-time representation of a `demes.Graph`

    :param yaml: The string representation of a demes model.
    :type yaml: str
    :param graph: The Python demes graph.
    :type graph: demes.Graph
    :param burnin: The number of generations of evolution to occur before the
                   events in the `yaml` take place.
    :type burnin: int

    The recommended method for construction is :meth:`fwdpy11.ForwardDemesGraph.from_demes`.

    .. versionadded:: 0.20.0
    """
    yaml: str = attr.ib()
    graph: demes.Graph = attr.ib()
    burnin: int = attr.ib()

    def __attrs_post_init__(self):
        super(ForwardDemesGraph, self).__init__(self.yaml, self.burnin)

    def _validate_pulses(graph: demes.Graph):
        from fwdpy11 import AmbiguousPulses

        unique_pulse_times = set([np.rint(p.time) for p in graph.pulses])
        for time in unique_pulse_times:
            pulses = [p for p in graph.pulses if np.rint(p.time) == time]
            dests = set()
            for p in pulses:
                if p.dest in dests:
                    raise AmbiguousPulses(
                        f"multiple pulse events into deme {p.dest} at time {time}"
                    )
                dests.add(p.dest)

    def number_of_demes(self) -> int:
        return len(self.graph.demes)

    @classmethod
    def from_demes(cls, model: typing.Union[str, demes.Graph], burnin: int):
        """
        Build from a demes graph.

        :param model: A `demes` model.
        :type model: str or demes.Graph
        :param burnin: The number of generations of evolution to occur before
                       the events in the `yaml` take place.
                       The value will be `burnin` times the total ancestral
                       population size.
        :type burnin: int
        """

        from ._functions.import_demes import _get_ancestral_population_size
        if isinstance(model, str) is False:
            yaml = str(model)
            graph = model
        else:
            yaml = model
            graph = demes.loads(yaml)
        cls._validate_pulses(graph)
        Nref = _get_ancestral_population_size(graph)
        burnin_time = int(np.rint(burnin*Nref))
        return cls(yaml, graph, burnin_time)

    @classmethod
    def tubes(cls, sizes: typing.List[int], burnin: int, *,
              burnin_is_exact: typing.Optional[bool] = None):
        total_popsize = sum(sizes)
        builder = demes.Builder()
        for i, j in enumerate(sizes):
            builder.add_deme(f"deme{i}", epochs=[dict(start_size=j)])
        graph = builder.resolve()
        if burnin_is_exact is not None and burnin is True:
            burnin_length = burnin
        else:
            burnin_length = burnin * total_popsize
        return cls.from_demes(graph, burnin_length)


def from_demes(
    dg: typing.Union[str, demes.Graph], burnin: int = 10
) -> "DemographicModelDetails":
    """
    Build a :class:`fwdpy11.demographic_models.DemographicModelDetails` object using demes.
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

    :rtype: :class:`fwdpy11.demographic_models.DemographicModelDetails`

    .. versionadded:: 0.14.0
    """
    from ._functions import demography_from_demes

    demog = demography_from_demes(dg, burnin)
    return demog
