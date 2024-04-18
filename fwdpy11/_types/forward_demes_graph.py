from dataclasses import dataclass
import decimal
import typing

import attr
import demes
import numpy as np
from deprecated import deprecated

from ..class_decorators import (
    attr_class_pickle_with_super,
    attr_class_to_from_dict,
)

import fwdpy11


def _round_via_decimal(value):
    with decimal.localcontext() as ctx:
        ctx.rounding = decimal.ROUND_HALF_UP
        return int(decimal.Decimal(value).to_integral_value())


@dataclass
class ForwardTimeInterval:
    """
    A half-open interval [start_time, end_time).
    The times represent values moving forward in time,
    in generations since the begining of a model.
    """

    start_time: int
    end_time: int


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
    :param round_non_integer_sizes: If `True`, a graph containing non-integer
                            values for epoch start and/or end sizes will
                            have those values rounded to the nearest integer.
    :type round_non_integer_sizes: bool

    The recommended method for construction is
    :meth:`ForwardDemesGraph.from_demes` or
    :meth:`ForwardDemesGraph.tubes`.

    .. versionadded:: 0.20.0
    """

    yaml: str = attr.ib()
    burnin: int = attr.ib()
    burnin_is_exact: int = attr.ib()
    round_non_integer_sizes: bool = attr.ib()
    graph: demes.Graph
    burnin_generation: int

    def __attrs_post_init__(self):
        self.graph = demes.loads(self.yaml)
        Nref = self._get_ancestral_population_size(self.graph)
        assert np.modf(Nref)[0] == 0.0
        if self.burnin_is_exact is True:
            burnin = self.burnin
        else:
            burnin = self.burnin * Nref
        self.burnin_generation = burnin
        super(ForwardDemesGraph, self).__init__(
            self.yaml, burnin, self.round_non_integer_sizes
        )
        x = self._sum_deme_sizes_at_time_zero()
        assert Nref == x, f"{Nref}, {x}, {self.yaml}"

    def __get_most_ancient_deme_start_time(self, dg: demes.Graph) -> demes.demes.Time:
        return max([d.start_time for d in dg.demes])

    def _get_ancestral_population_size(self, dg: demes.Graph) -> int:
        """
        Need this for the burnin time.

        If there are > 1 demes with the same most ancient start_time,
        then the ancestral size is considered to be the size
        of all those demes (size of ancestral metapopulation).
        """
        oldest_deme_time = self.__get_most_ancient_deme_start_time(dg)

        rv = sum(
            [
                _round_via_decimal(e.start_size)
                for d in dg.demes
                for e in d.epochs
                if e.start_time == oldest_deme_time
            ]
        )
        if rv == 0:
            raise ValueError("could not determinine ancestral metapopulation size")
        return rv

    def _minimal_end_time_from_demes_graph(self):
        return min(
            [_round_via_decimal(i.end_time) for i in self.graph.in_generations().demes]
        )

    def number_of_demes(self) -> int:
        return len(self.graph.demes)

    @classmethod
    def from_demes(
        cls,
        model: typing.Union[str, demes.Graph],
        burnin: int,
        *,
        burnin_is_exact: typing.Optional[bool] = None,
        round_non_integer_sizes: typing.Optional[bool] = None,
    ):
        """
        Build from a demes graph.

        :param model: A `demes` model.
        :type model: str or demes.Graph
        :param burnin: The number of generations of evolution to occur before
                       the events in the `yaml` take place.
        :type burnin: int
        :param burnin_is_exact: If `False`, `burnin` will be treated as an
                                integer multiple of the sum of ancestor
                                population sizes. If `True`, `burnin` will
                                be used as-is.
        :type burnin_is_exact: bool
        :type round_non_integer_sizes: If `True`, a graph containing
                                non-integer values for epoch start
                                and/or end sizes will have those values
                                rounded to the nearest integer.

        :returns: A forward-time representation of a demes graph
        :rtype: fwdpy11.ForwardDemesGraph
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
            round_sizes = False
        else:
            round_sizes = round_non_integer_sizes

        return cls(yaml, burnin, exact, round_sizes)

    @classmethod
    def tubes(
        cls,
        sizes: typing.List[int],
        burnin: int,
        *,
        burnin_is_exact: typing.Optional[bool] = None,
        round_non_integer_sizes: typing.Optional[bool] = None,
    ):
        """
        Build a demographic model for one or more demes of constant size.

        :param sizes: Deme sizes
        :type sizes: typing.List[int]
        :param burnin: The number of generations of evolution to simulate.
        :type burnin: int
        :param burnin_is_exact: If `False`, `burnin` will be treated as an
                                integer multiple of the sum of ancestor
                                population sizes. If `True`, `burnin` will
                                be used as-is.
        :type burnin_is_exact: bool
        :type round_non_integer_sizes: If `True`, a graph containing
                                non-integer values for epoch start
                                and/or end sizes will have those values
                                rounded to the nearest integer.

        :returns: A forward-time representation of a demes graph
        :rtype: fwdpy11.ForwardDemesGraph
        """
        builder = demes.Builder()
        for i, j in enumerate(sizes):
            builder.add_deme(f"deme{i}", epochs=[dict(start_size=j)])
        graph = builder.resolve()
        return cls.from_demes(
            graph,
            burnin,
            burnin_is_exact=burnin_is_exact,
            round_non_integer_sizes=round_non_integer_sizes,
        )

    def to_backwards_time(self, forwards_time: int) -> typing.Optional[int]:
        """
        Convert an integer value representing time moving in the
        forward direction to time in the past ("backwards time").

        If no conversion is possible, `None` is returned.

        :param forwards_time: The value to convert
        :type forwards_time: int

        :returns: The converted time
        :type: typing.Optional[int]
        """
        if forwards_time >= 0 and forwards_time <= self.final_generation:
            minimal_end_time = self._minimal_end_time_from_demes_graph()
            return self.final_generation - forwards_time + minimal_end_time
        return None

    def to_forwards_time(self, backwards_time: int) -> typing.Optional[int]:
        """
        Convert an integer value representing time moving in the
        backward direction to time in the forward direction.

        If no conversion is possible, `None` is returned.

        :param backwards_time: The value to convert
        :type backwards_time: int

        :returns: The converted time
        :type: typing.Optional[int]
        """
        if backwards_time >= 0 and self.final_generation - backwards_time >= 0:
            minimal_end_time = self._minimal_end_time_from_demes_graph()
            return self.final_generation - backwards_time + minimal_end_time
        return None

    def deme_time_intervals(self) -> typing.List[ForwardTimeInterval]:
        """
        :returns: a list of the time intervals during which each deme exists.
        :rtype: list[fwdpy11.ForwardTimeInterval]

        The order of the demes is the same as in the input :class:`demes.Graph`.
        """
        rv = []
        graph_in_generations = self.graph.in_generations()
        for deme in graph_in_generations.demes:
            start = deme.start_time
            if start == float("inf"):
                start = 0
            else:
                # The +1 takes the 1/2 open backwards time
                # interval and gives us the first generation
                # when individuals are born in this deme
                start = self.to_forwards_time(_round_via_decimal(start)) + 1
            # The +1 makes it a 1/2-open interval
            end = self.to_forwards_time(_round_via_decimal(deme.end_time)) + 1
            rv.append(ForwardTimeInterval(start, end))
        return rv

    def epoch_time_intervals(self) -> typing.List[typing.List[ForwardTimeInterval]]:
        """
        :returns: a list of the time intervals during which
                  each epoch of each deme exists.
        :rtype: list[list[fwdpy11.ForwardTimeInterval]]

        The order of the demes is the same as in the input :class:`demes.Graph`.
        Within each deme, the epoch order is the same as in the input graph.
        """
        rv = []
        graph_in_generations = self.graph.in_generations()
        for deme in graph_in_generations.demes:
            temp = []
            for epoch in deme.epochs:
                start = epoch.start_time
                if start == float("inf"):
                    start = 0
                else:
                    # The +1 takes the 1/2 open backwards time
                    # interval and gives us the first generation
                    # when individuals are born in this deme
                    start = self.to_forwards_time(_round_via_decimal(start)) + 1
                # The +1 makes it a 1/2-open interval
                end = self.to_forwards_time(_round_via_decimal(epoch.end_time)) + 1
                temp.append(ForwardTimeInterval(start, end))
            rv.append(temp)
        return rv

    @property
    def initial_sizes(self) -> typing.List:
        """
        A list of the nonzero (parental) deme sizes at time 0 of the model
        (thinking forwards in time).

        Due to the structure of a demes graph, the indexes of the list
        are the deme integer identifiers.

        This property is useful for setting up a
        :class:`fwdpy11.DiploidPopulation`.
        """
        return [i for i in self._parental_deme_sizes_at_time_zero() if i > 0]

    @property
    @deprecated(reason="prefer initial_sizes")
    def initial_sizes_list(self) -> typing.List:
        return self.initial_sizes

    @property
    def final_generation(self) -> int:
        """
        This is the final offspring generation of the model.
        In other words, a model evolved from parental generation
        zero until time zero of the graph (backwards in time)
        will result in individuals whose birth times equal
        the value of this property.
        """
        return self._model_end_time() - 1

    @property
    def burnin_duration(self) -> int:
        """
        The length of the burning time.
        Equal to `burnin_generation` + 1.
        """
        return self.burnin_generation + 1

    @property
    @deprecated(reason="prefer final_generation")
    def total_simulation_length(self) -> int:
        return self.final_generation

    @property
    def model_duration(self) -> int:
        """
        The duration of the model from time 0 until the
        most recent time point.
        """
        return self.final_generation

    @property
    def deme_labels(self) -> typing.Dict[int, str]:
        """
        A dictionary mapping integer deme IDs to deme names
        (strings).
        """
        return {i: j.name for i, j in enumerate(self.graph.demes)}

    @property
    def demes_at_final_generation(self) -> typing.List[int]:
        """
        List of integer ids of extant demes in the final generation
        of the demes graph.
        """
        epoch_end_times = []
        for deme in self.graph.demes:
            for epoch in deme.epochs:
                epoch_end_times.append(epoch.end_time)
        min_end_time = min(epoch_end_times)
        rv = []
        for deme in self.deme_labels.keys():
            for epoch in self.graph.demes[deme].epochs:
                if epoch.end_time == min_end_time:
                    rv.append(deme)
        return rv

    @property
    def demes_graph(self) -> demes.Graph:
        """
        Obtain the internal representation of the :class:`demes.Graph`

        This object is not necessarily equivalent to the input graph:

        * It is fully resolved.
        * It has been rounded to discrete time and deme sizes.
        """
        return demes.loads(self._demes_graph())
