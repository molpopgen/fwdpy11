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
    :param burnin: The number of generations of evolution to occur before the
                   events in the `yaml` take place.
    :type burnin: int
    :param burnin_is_exact: If `False`, `burnin` will be treated as an integer
                            multiple of the sum of ancestor population sizes.
                            If `True`, `burnin` will be used as-is.
    :type burnin_is_exact: bool
    :type round_non_integer_sizes: If `True`, a graph containing non-integer
                            values for epoch start and/or end sizes will
                            have those values rounded to the nearest integer.

    The recommended method for construction is
    :meth:`fwdpy11.ForwardDemesGraph.from_demes` or
    :meth:`fwdpy11.ForwardDemesGraph.tubes`

    .. versionadded:: 0.20.0
    """
    yaml: str = attr.ib()
    burnin: int = attr.ib()
    burnin_is_exact: int = attr.ib()
    round_non_integer_sizes: bool = attr.ib()
    graph: demes.Graph

    def __attrs_post_init__(self):
        from ._functions.import_demes import _get_ancestral_population_size
        self.graph = demes.loads(self.yaml)
        if self.round_non_integer_sizes is True:
            self.graph = ForwardDemesGraph._round_deme_sizes(self.graph)

        ForwardDemesGraph._validate_pulses(self.graph)
        updated_yaml = str(self.graph)
        Nref = _get_ancestral_population_size(self.graph)
        if self.burnin_is_exact is True:
            burnin = self.burnin
        else:
            burnin = int(np.rint(self.burnin)*Nref)
        super(ForwardDemesGraph, self).__init__(updated_yaml, burnin)

    def _validate_pulses(graph: demes.Graph):
        unique_pulse_times = set([np.rint(p.time) for p in graph.pulses])
        for time in unique_pulse_times:
            pulses = [p for p in graph.pulses if np.rint(p.time) == time]
            dests = set()
            for p in pulses:
                if p.dest in dests:
                    warnings.warn(
                        f"multiple pulse events into deme {p.dest} at time {time}."
                        + " The effect of these pulses will depend on the order "
                        + "in which they are applied."
                        + "To avoid unexpected behavior, "
                        + "the graph can instead be structured to"
                        + " introduce a new deme at this time with"
                        + " the desired ancestry proportions or to specify"
                        + " concurrent pulses with multiple sources.",
                        UserWarning)
                dests.add(p.dest)

    def _reject_non_integer_sizes(graph: demes.Graph):
        for deme in graph.demes:
            for i, epoch in enumerate(deme.epochs):
                for size in [epoch.start_size, epoch.end_size]:
                    if np.isfinite(size) and np.modf(size)[0] != 0.0:
                        raise ValueError(
                            f"deme {deme.name} has non-integer size {size} in epoch {i}")

    def _round_deme_sizes(graph: demes.Graph) -> demes.Graph:
        import copy
        graph_copy = copy.deepcopy(graph)
        for deme in graph_copy.demes:
            for epoch in deme.epochs:
                epoch.start_size = int(np.rint(epoch.start_size))
                epoch.end_size = int(np.rint(epoch.end_size))
                for i in [epoch.start_size, epoch.end_size]:
                    if i <= 0.0:
                        raise ValueError(
                            "rounding resulted in a deme size of 0")

        return graph_copy

    def number_of_demes(self) -> int:
        return len(self.graph.demes)

    @classmethod
    def from_demes(cls, model: typing.Union[str, demes.Graph], burnin: int, *,
                   burnin_is_exact: typing.Optional[bool] = None,
                   round_non_integer_sizes: typing.Optional[bool] = None):
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
        # NOTE: the logic is messy and we are losing
        # track of "provenance" of the demes input
        if isinstance(model, str) is False:
            yaml = str(model)
        else:
            yaml = model

        if burnin_is_exact is None:
            exact = False
        else:
            exact = burnin_is_exact

        if round_non_integer_sizes is None:
            round_sizes = True
        else:
            round_sizes = round_non_integer_sizes

        return cls(yaml, burnin, exact, round_sizes)

    @classmethod
    def tubes(cls, sizes: typing.List[int], burnin: int, *,
              burnin_is_exact: typing.Optional[bool] = None,
              round_non_integer_sizes: typing.Optional[bool] = None):
        builder = demes.Builder()
        for i, j in enumerate(sizes):
            builder.add_deme(f"deme{i}", epochs=[dict(start_size=j)])
        graph = builder.resolve()
        return cls.from_demes(graph, burnin,
                              burnin_is_exact=burnin_is_exact,
                              round_non_integer_sizes=round_non_integer_sizes)


def from_demes(
    dg: typing.Union[str, demes.Graph], burnin: int = 10, *,
    round_non_integer_sizes: typing.Optional[bool] = None
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

    demog = demography_from_demes(
        dg, burnin, round_non_integer_sizes=round_non_integer_sizes)
    return demog
