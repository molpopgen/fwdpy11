#
# Copyright (C) 2017 Kevin Thornton <krthornt@uci.edu>
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

import abc
import typing
import warnings

import attr
import fwdpy11
import numpy as np
from fwdpy11.class_decorators import (
    attr_class_to_from_dict_no_recurse,
)

from fwdpy11._types.demographic_model_details import DemographicModelDetails
from fwdpy11._types.forward_demes_graph import ForwardDemesGraph


@attr_class_to_from_dict_no_recurse
@attr.s(kw_only=True, frozen=True, slots=True)
class MutationAndRecombinationRates(object):
    """
    Stores and validates the mutation and recombination rate parameters
    of a simulation.

    Instances of this class are created by ``kwargs`` that populate
    attributes of the same name:

    :param neutral_mutation_rate:
    :type neutral_mutation_rate: float
    :param selected_mutation_rate:
    :type selected_mutation_rate: float
    :param recombination_rate:
    :type recombination_rate: float or None

    Instances of this class are passed as the ``rates`` ``kwarg``
    to :class:`fwdpy11.ModelParams`.

    .. versionadded:: 0.8.0
    """

    neutral_mutation_rate: float = attr.ib(converter=float)
    selected_mutation_rate: float = attr.ib(converter=float)
    recombination_rate = attr.ib()

    @neutral_mutation_rate.validator
    @selected_mutation_rate.validator
    def validate_individual_rate(self, attribute, value):
        if not np.isfinite(value):
            raise ValueError(f"{attribute} must be finite, but we got {value} instead")
        if value < 0.0:
            raise ValueError(f"{attribute} must be >= 0.0, but we got {value} instead")

    @recombination_rate.validator
    def validate_recombination_rate(self, attribute, value):
        if value is None:
            return
        try:
            value = float(value)
        except TypeError:
            raise ValueError(
                f"{attribute} must be convertible to float, but we got {value} instead"
            )
        attr.validators.instance_of(float)(self, attribute, value)
        self.validate_individual_rate(attribute, value)


def _convert_rates(value):
    if isinstance(value, MutationAndRecombinationRates):
        return value

    try:
        return MutationAndRecombinationRates(
            neutral_mutation_rate=value[0],
            selected_mutation_rate=value[1],
            recombination_rate=value[2],
        )
    except KeyError:
        return MutationAndRecombinationRates(**value)


@typing.runtime_checkable
class ValidateEventTimings(typing.Protocol):
    @abc.abstractmethod
    def validate_timings(self, deme: int, demography: ForwardDemesGraph) -> None:
        raise NotImplementedError


def validate_timings(
    gvalue: ValidateEventTimings, deme: int, demography: ForwardDemesGraph
) -> None:
    assert isinstance(gvalue, ValidateEventTimings)
    return gvalue.validate_timings(deme, demography)


@attr_class_to_from_dict_no_recurse
@attr.s(kw_only=True, frozen=True, slots=True)
class ModelParams(object):
    """
    This class stores and validates the parameters of a simulation.
    Instances of this class are constructed using ``kwargs``
    and instance attributes are immutable once initialized.

    This class accepts the following ``kwargs``, which are also
    the attribute names for instances:

    :param nregions: List of regions where neutral mutations occur
    :type nregions: list[fwdpy11.Region]
    :param sregions: List of regions where selected mutations occur
    :type sregions: list[fwdpy11.Sregion]
    :param recregions: List of regions where recombination events occur
    :type recregions: list[fwdpy11.Region] or list[object]
    :param rates: The neutral mutation rate, selected mutation rate, and
                  total recombination rate, respectively.
                  See below for more details.
    :type rates: list or fwdpy11.MutationAndRecombinationRates
    :param demography: The demographic model to simulate
    :type demography: object
    :param simlen: The number of time steps to evolve
    :type simlen: int
    :param prune_selected: If ``True``, remove selected fixations from
                           the population when they are first detected
                           as fixed.
    :type prune_selected: bool
    :param allow_residual_selfing: If ``True``, then an individual may
                                   be chosen twice as a parent.
                                   If ``False``, two distinct parents
                                   are required.
    :type allow_residual_selfing: bool

    .. warning::

        When ``allow_residual_selfing`` is ``False``, an exception
        will be raised if the parental deme has a size of 1.
        The reason is that outcrossing is not possible in order
        to generate the offspring.

    .. note::

        To initialize the ``rates`` field, we require an instance of
        fwdpy11.MutationAndRecombinationRates or a list of length three (3)
        specifying the three rates.
        The two mutation rates must be non-negative floats. For
        the recombination rate, the third value must also
        be a non-negative float if all objects in ``recrates``
        are instances of :class:`fwdpy11.Region`. However,
        if they are instead instances of any mixture of
        :class:`fwdpy11.PoissonInterval`, :class:`fwdpy11.PoissonPoint`,
        :class:`fwdpy11.BinomialPoint`, or :class:`fwdpy11.BinomialInterval`,
        then the final element in ``rates`` must be ``None``.
        See the :ref:`section <geneticmaps_vignette>` on setting
        recombination rates for details.


    If ``prune_selected`` is ``True``, fixations will be remove
    during the simulation. This removal is a sensible thing to do
    if individual relative fitnesses are constant up to a multiplicative
    constant.

    In a multiplicative model, removing a fixation would be equivalent to dividing all values
    by the same value, which does not affect their *relative*
    values.

    However, removing a fixation in an additive model *does* change
    the relative values.
    Likewise, mappings from genetic value to fitness such as Gaussian
    stabilizing selection, or adding random noise to genetic values,
    also changes relative values.

    This class will raise a warning if ``prune_selected`` is ``True``
    and there is a possibility that removing fixations will change
    relative fitness.
    Such warnings should not be ignored as they suggest that a simulation
    may not be giving biologically sensible results.

    .. versionadded:: 0.1.1

    .. versionchanged:: 0.2.0
        Changed this from a horrible class hierarchy
        into a much simpler, single class.

    .. versionchanged:: 0.6.0

        Updated to support `fwdpy11.DiscreteDemography`

    .. versionchanged:: 0.8.0

        Refactored class internals using ``attrs``.
        Mutation and recombination rates now
        stored in :class:`fwdpy11.MutationAndRecombinationRates`

    .. versionchanged:: 0.14.0

        Remove deprecated kwargs pself and popsizes.

    .. versionchanged:: 0.20.0

        Add `allow_residual_selfing`
    """

    nregions = attr.ib(factory=list)
    sregions = attr.ib(factory=list)
    recregions = attr.ib(factory=list)
    rates: MutationAndRecombinationRates = attr.ib(converter=_convert_rates)
    gvalue = attr.ib(default=None)
    demography: typing.Union[ForwardDemesGraph, DemographicModelDetails] = attr.ib(
        default=None
    )
    simlen: int = attr.ib(converter=int, default=0)
    prune_selected: bool = attr.ib(default=True)
    allow_residual_selfing: bool = attr.ib(default=True)

    @nregions.validator
    def validate_nregions(self, attribute, value):
        for i in value:
            attr.validators.instance_of(fwdpy11.Region)(self, attribute, i)

    @sregions.validator
    def validate_sregions(self, attribute, value):
        for i in value:
            attr.validators.instance_of(fwdpy11.Sregion)(self, attribute, i)
            try:
                if i.shape != self.gvalue.shape:
                    e = "Sregion and genetic value "
                    "dimension mismatch: {} {}, {} {}".format(
                        type(i), i.shape, type(self.gvalue), self.gvalue.shape
                    )
                    raise ValueError(e)
            except AttributeError:
                for g in self.gvalue:
                    if i.shape != g.shape:
                        e = "Sregion and genetic value "
                        "dimension mismatch: {} {}, {} {}".format(
                            type(i), i.shape, type(self.gvalue), g.shape
                        )
                        raise ValueError(e)

    @recregions.validator
    def validate_recregions(self, attribute, value):
        if len(value) > 0 and all([isinstance(i, fwdpy11.Region) for i in value]):
            warnings.warn(
                "using a list of Regions for recregions is deprecated",
                UserWarning,
                stacklevel=2,
            )
        try:
            for i in value:
                attr.validators.instance_of(fwdpy11.Region)(self, attribute, i)
        except TypeError:
            try:
                for i in value:
                    valid = False
                    for k in [
                        fwdpy11.PoissonCrossoverGenerator,
                        fwdpy11.NonPoissonCrossoverGenerator,
                    ]:
                        if isinstance(i, k):
                            valid = True
                            break

                    if not valid:
                        raise TypeError(f"invalid recregion type {type(i)}")
                if not all([i.discrete for i in value]) and not all(
                    [not i.discrete for i in value]
                ):
                    warnings.warn(
                        "genetic map has a mix of discrete=True" "and discrete=False",
                        stacklevel=2,
                    )

            except TypeError:
                raise

    @rates.validator
    def rates_validator(self, attribute, value):
        if value.recombination_rate is None:
            for i in self.recregions:
                if not isinstance(
                    i, fwdpy11.PoissonCrossoverGenerator
                ) and not isinstance(i, fwdpy11.NonPoissonCrossoverGenerator):
                    raise ValueError(
                        f"recombination rate of {value.recombination_rate}"
                        " must be paired with"
                        " instances of fwdpy11.GeneticMapUnit"
                    )
        else:
            for i in self.recregions:
                if not isinstance(i, fwdpy11.Region):
                    raise ValueError(
                        f"recombination rate of {value.recombination_rate}"
                        " must be paired with"
                        " instances of fwdpy11.Region"
                    )

    @gvalue.validator
    def validate_gvalue(self, attribute, value):
        try:
            for deme, i in enumerate(value):
                attr.validators.instance_of(fwdpy11.DiploidGeneticValue)(
                    self, attribute, i
                )
                validate_timings(i, deme, self.demography)
        except TypeError:
            try:
                attr.validators.instance_of(fwdpy11.DiploidGeneticValue)(
                    self, attribute, value
                )
                validate_timings(value, 0, self.demography)
            except TypeError:
                raise

    @demography.validator
    def validate_demography(self, attribute, value):
        if value is None:
            raise ValueError("no demographic model specified")

        if isinstance(value, fwdpy11.ForwardDemesGraph):
            return

        if isinstance(value, DemographicModelDetails):
            return
        raise ValueError(f"Unknown type for {attribute}: {type(value)}")

    @simlen.validator
    def validate_simlen(self, attribute, value):
        if value <= 0:
            raise ValueError(f"{attribute} must be >= 0, but we got {value} instead")

    def _is_multiplicative(self, gvalue):
        multiplicative = True

        if not isinstance(gvalue, fwdpy11.Multiplicative):
            multiplicative = False
        if multiplicative is True:
            if gvalue.gvalue_to_fitness is not None:
                if not isinstance(
                    gvalue.gvalue_to_fitness, fwdpy11.GeneticValueIsFitness
                ):
                    multiplicative = False
            if gvalue.noise is not None and not isinstance(
                gvalue.noise, fwdpy11.NoNoise
            ):
                multiplicative = False

        if not multiplicative and self.prune_selected is True:
            msg = f"""
            the genetic value {gvalue} may not be proportional to fitness
            up to a multiplicative constant. Because prune_selected=True,
            the removal of fixations during simulation may lead to unexpected
            or incorrect results. See the ModelParams documentation
            for discussion.
            """
            warnings.warn(fwdpy11.PruneSelectedWarning(msg))

    def __attrs_post_init__(self):
        if isinstance(self.gvalue, list):
            for g in self.gvalue:
                self._is_multiplicative(g)
        else:
            self._is_multiplicative(self.gvalue)
