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
    NEVER = 0
    DURATION = 1
    COMPLETION = 2


def non_negative_value(_, attribute, value):
    if value < 0:
        raise ValueError(f"{attribute.name} must be >= 0, got {value}")


def greater_than_zero(_, attribute, value):
    if value < 0:
        raise ValueError(f"{attribute.name} must be > 0")


def maximum_greater_minimum(instance, attribute, value):
    if value <= instance.minimum:
        raise ValueError(f"{attribute.name} must be > minimum value")


def right_greater_left(instance, attribute, value):
    if value <= instance.left:
        raise ValueError(f"{attribute.name} must be > left value")


def is_polymorphic_frequency(_, attribute, value):
    if not 0.0 < value < 1.0:
        raise ValueError(f"{attribute.name} must be 0.0 < {attribute.name} < 1.0")


@attr.s(auto_attribs=True, frozen=True)
class AlleleCount:
    count: int = attr.ib(
        validator=[attr.validators.instance_of(int), greater_than_zero]
    )


@attr.s(auto_attribs=True, frozen=True)
class AlleleCountRange:
    minimum: int = attr.ib(
        validator=[attr.validators.instance_of(int), greater_than_zero]
    )
    maximum: int = attr.ib(
        validator=[
            attr.validators.instance_of(int),
            greater_than_zero,
            maximum_greater_minimum,
        ]
    )


@attr.s(auto_attribs=True, frozen=True)
class FrequencyRange:
    minimum: float = attr.ib(
        validator=[attr.validators.instance_of(float), is_polymorphic_frequency]
    )
    maximum: float = attr.ib(
        validator=[
            attr.validators.instance_of(float),
            is_polymorphic_frequency,
            maximum_greater_minimum,
        ]
    )


@attr.s(auto_attribs=True, kw_only=True, frozen=True)
class PositionRange:
    left: float = attr.ib(
        validator=[attr.validators.instance_of(float), non_negative_value]
    )
    right: float = attr.ib(
        validator=[
            attr.validators.instance_of(float),
            greater_than_zero,
            right_greater_left,
        ]
    )


@attr.s(auto_attribs=True, kw_only=True, frozen=True)
class NewMutationParameters:
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
    should_terminate: bool = attr.ib(validator=attr.validators.instance_of(bool))
    condition_met: bool = attr.ib(validator=attr.validators.instance_of(bool))


@attr.s(auto_attribs=True, kw_only=True)
class ConditionalModelOutput:
    pop: fwdpy11.DiploidPopulation
    index: int
    num_nodes: int


class GlobalFixation(object):
    def __call__(
        self, pop: fwdpy11.DiploidPopulation, index: int, key: tuple
    ) -> SimulationStatus:
        if pop.mutations[index].key != key:
            return SimulationStatus(False, False)
        if pop.mcounts[index] == 0:
            return SimulationStatus(True, False)
        if pop.mcounts[index] == 2 * pop.N:
            return SimulationStatus(False, True)
        return SimulationStatus(False, False)


@attr.s(auto_attribs=True, frozen=True)
class FocalDemeFixation:
    deme: int = attr.ib(
        validator=[attr.validators.instance_of(int), non_negative_value]
    )

    def __call__(self, pop: fwdpy11.DiploidPopulation, index, key) -> SimulationStatus:
        deme_sizes = pop.deme_sizes(as_dict=True)
        if pop.mutations[index].key != key:
            return SimulationStatus(False, False)
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
    simplification_interval: int = attr.ib(
        validator=[attr.validators.instance_of(int), non_negative_value], default=100
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
