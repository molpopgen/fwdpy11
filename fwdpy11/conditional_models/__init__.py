#
# Copyright (C) 2021 Kevin Thornton <krthornt@uci.edu>
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

"""
Defining models with specific outcomes.

The API of this module may be subjet to change.
Hence, the documentation is minimal.
"""

import typing
from enum import Enum

import attr
import fwdpy11
from fwdpy11.class_decorators import attr_class_to_from_dict


class AddMutationFailure(Exception):
    pass


class OutOfAttempts(Exception):
    pass


class AncientSamplePolicy(Enum):
    """
    When to record "preserved" or "ancient" samples:
    """

    NEVER = 0
    """
    Never.
    """
    DURATION = 1
    """
    During the entire sojourn of the mutation from when it is first
    added to when the terminiation condition is first met.
    """
    COMPLETION = 2
    """
    Only record the generation when the terminiation condition is first met.
    """


def _non_negative_value(_, attribute, value):
    if value < 0:
        raise ValueError(f"{attribute.name} must be >= 0, got {value}")


def _greater_than_zero(_, attribute, value):
    if value < 0:
        raise ValueError(f"{attribute.name} must be > 0")


def _maximum_greater_minimum(instance, attribute, value):
    if value <= instance.minimum:
        raise ValueError(f"{attribute.name} must be > minimum value")


def _right_greater_left(instance, attribute, value):
    if value <= instance.left:
        raise ValueError(f"{attribute.name} must be > left value")


def _is_polymorphic_frequency(_, attribute, value):
    if not 0.0 < value < 1.0:
        raise ValueError(
            f"{attribute.name} must be 0.0 < {attribute.name} < 1.0")


@attr.s(auto_attribs=True, frozen=True)
class AlleleCount:
    """
    Specify a *number* of copies of a mutation.

    :param count: Initial number of copies of a mutation.
                  This value must be `> 0`.
    :type count: int
    """

    count: int = attr.ib(
        validator=[attr.validators.instance_of(int), _greater_than_zero]
    )


@attr.s(auto_attribs=True, frozen=True)
class AlleleCountRange:
    """
    Specify a range for a *number* of copies of a mutation.

    :param minimum: Minimum number of copies of a mutation.
                    This value must be `> 0`.
    :type minimum: int
    :param maximum: Maximum number of copies of a mutation.
                    This value must be `> minimum`.
    :type maximum: int
    """

    minimum: int = attr.ib(
        validator=[attr.validators.instance_of(int), _greater_than_zero]
    )
    maximum: int = attr.ib(
        validator=[
            attr.validators.instance_of(int),
            _greater_than_zero,
            _maximum_greater_minimum,
        ]
    )


@attr.s(auto_attribs=True, frozen=True)
class FrequencyRange:
    """
    Specify a range for the initial *frequency* of a mutation

    :param minimum: Minimum frequency of a mutation.
                    This value must be `> 0.0` and `< 1.0`.
    :type minimum: float
    :param maximum: Maximum frequency of a mutation.
                    This value must be `> minimum`.
    :type maximum: float
    """

    minimum: float = attr.ib(
        validator=[attr.validators.instance_of(
            float), _is_polymorphic_frequency]
    )
    maximum: float = attr.ib(
        validator=[
            attr.validators.instance_of(float),
            _is_polymorphic_frequency,
            _maximum_greater_minimum,
        ]
    )


@attr.s(auto_attribs=True, kw_only=True, frozen=True)
class PositionRange:
    """
    Specify a half-open interval within which to place a mutation.

    :param left: The left edge, inclusive, of the interval
    :type left: float
    :param right: The right edge, exclusive, of the interval
    :type right: float
    """

    left: float = attr.ib(
        validator=[attr.validators.instance_of(float), _non_negative_value]
    )
    right: float = attr.ib(
        validator=[
            attr.validators.instance_of(float),
            _greater_than_zero,
            _right_greater_left,
        ]
    )


@attr.s(auto_attribs=True, kw_only=True, frozen=True)
class NewMutationParameters:
    """
    Details of a new mutation to add to a population.

    Class instances are created via keyword arguments that
    become attribute names:

    :param deme: The id of the deme in which to add the mutation
    :type deme: Optional[int]
    :param frequency: The frequency of the new mutation.
    :type frequency: Union[AlleleCount, AlleleCountRange, FrequencyRange]
    :param position: Where to put the new mutation
    :type position: PositionRange
    :param data: The specifics of the new mutation
    :type data: :class:`fwdpy11.NewMutationData`
    """

    deme: typing.Optional[int] = None
    frequency: typing.Union[AlleleCount, AlleleCountRange, FrequencyRange] = attr.ib(
        validator=attr.validators.instance_of(
            (AlleleCount, AlleleCountRange, FrequencyRange)
        )
    )
    position: PositionRange = attr.ib(
        validator=attr.validators.instance_of(PositionRange)
    )
    data: fwdpy11.NewMutationData = attr.ib(
        validator=attr.validators.instance_of(fwdpy11.NewMutationData)
    )


@attr.s(frozen=True)
class SimulationStatus:
    """
    The return value of a stopping condition callable.

    :param should_terminate: Set to `True` if the simulation should be terminated.
    :param condition_met: Set to `True` if the stopping condition has been met.

    For examples, see implementations of :class:`GlobalFixation` and :class:`FocalDemeFixation`.
    """

    should_terminate: bool = attr.ib(
        validator=attr.validators.instance_of(bool))
    condition_met: bool = attr.ib(validator=attr.validators.instance_of(bool))


@attr.s(auto_attribs=True, kw_only=True)
class ConditionalModelOutput:
    pop: fwdpy11.DiploidPopulation
    """The population with the added mutation"""
    params: fwdpy11.ModelParams
    """The evolved model parameters"""
    mutation_index: int
    """The index of the new mutation in pop.mutations"""
    num_descendant_nodes: int
    """The number of alive nodes initially containing the new mutation"""


class GlobalFixation(object):
    """
    A terminiation condition monitoring for global fixation of
    the mutation.
    """

    def __call__(
        self, pop, index: int, key: tuple
    ) -> SimulationStatus:
        if pop.mutations[index].key != key:
            # The key has changed, meaning the mutation is
            # flagged for recycling.
            # First, check if it is in the fixations list
            for m in pop.fixations:
                if m.key == key:
                    # It is fixed, so we are done
                    return SimulationStatus(True, False)
            # The mutation is gone from the simulation
            return SimulationStatus(True, False)
        if pop.mcounts[index] == 0:
            return SimulationStatus(True, False)
        if pop.mcounts[index] == 2 * pop.N:
            return SimulationStatus(False, True)
        return SimulationStatus(False, False)


@attr.s(auto_attribs=True, frozen=True)
class FocalDemeFixation:
    """
    A terminiation condition checking for fixation in a specific deme.
    """

    deme: int = attr.ib(
        validator=[attr.validators.instance_of(int), _non_negative_value]
    )

    def __call__(self, pop, index, key) -> SimulationStatus:
        deme_sizes = pop.deme_sizes(as_dict=True)
        if pop.mutations[index].key != key:
            # check for a global fixation
            for m in pop.fixations:
                if m.key == key:
                    return SimulationStatus(True, False)
            # The mutation is gone from the simulation
            return SimulationStatus(True, False)
        if self.deme not in deme_sizes:
            return SimulationStatus(False, False)
        count = 0

        # TODO: we can probably do better here
        # by using the numpy interface to the metadata
        for md in pop.diploid_metadata:
            if md.deme == self.deme:
                for genome in [
                    pop.diploids[md.label].first,
                    pop.diploids[md.label].second,
                ]:
                    if index in pop.haploid_genomes[genome].smutations:
                        count += 1

        if count == 0:
            return SimulationStatus(True, False)
        if count == 2 * deme_sizes[self.deme]:
            return SimulationStatus(False, True)
        return SimulationStatus(False, False)


@attr_class_to_from_dict
@attr.s(frozen=True)
class EvolveOptions:
    """
    Options to pass on too :func:`fwdpy11.evolvets`.
    """

    simplification_interval: int = attr.ib(
        validator=[attr.validators.instance_of(int), _non_negative_value], default=100
    )
    suppress_table_indexing: bool = attr.ib(
        validator=attr.validators.instance_of(bool), default=True
    )
    record_gvalue_matrix: bool = attr.ib(
        validator=attr.validators.instance_of(bool), default=False
    )
    preserve_first_generation: bool = attr.ib(
        validator=attr.validators.instance_of(bool), default=False
    )


from ._selective_sweep import _selective_sweep as selective_sweep  # NOQA
from ._track_added_mutation import _track_added_mutation as track_added_mutation  # NOQA
