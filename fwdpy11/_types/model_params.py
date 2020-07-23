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
import warnings

import attr
import numpy as np

import fwdpy11
from fwdpy11.class_decorators import (attr_add_asblack,
                                      attr_class_to_from_dict_no_recurse)

_common_attr_attribs = {
    "kw_only": True,
    "frozen": True,
    "slots": True,
    "repr_ns": "fwdpy11",
}


@attr_add_asblack
@attr_class_to_from_dict_no_recurse
@attr.s(**_common_attr_attribs)
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


def _convert_demography(o):
    try:
        o.items()
        return fwdpy11.DiscreteDemography(**o)
    except AttributeError:
        return o


@attr_add_asblack
@attr_class_to_from_dict_no_recurse
@attr.s(**_common_attr_attribs)
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
    :type recregions: list[fwdpy11.Region] or list[fwdpy11.GeneticMapUnit]
    :param rates: The neutral mutation rate, selected mutation rate, and
                  total recombination rate, respectively.
                  See below for more details.
    :type rates: list or fwdpy11.MutationAndRecombinationRates
    :param demography: The demographic model to simulate
    :type demography: fwdpy11.DiscreteDemography
    :param simlen: The number of time steps to evolve
    :type simlen: int
    :param prune_selected: If ``True``, remove selected fixations from
                           the population when they are first detected
                           as fixed.
    :type prune_selected: bool

    The following ``kwargs``/attributes are pending deprecation along
    with simulations without tree sequence recording:

    :param pself: The probabilility that an individual selfs.
                  For simulations with tree sequence recording,
                  this parameter has no effect.  Instead, use
                  the methods described :ref:`here <softselection>`.
    :type self: float
    :param popsizes: A list of population sizes over time.
                     For simulations with tree sequence recording,
                     this parameter has no effect.  Instead, use
                     the methods described :ref:`here <softselection>`.
    :type popsizes: numpy.ndarray

    .. note::

        To initialize the ``rates`` field, we require an instance of
        fwdpy11.MutationAndRecombinationRates or a list of length three (3)
        specifying the three rates.
        The two mutation rates must be non-negative floats. For
        the recombination rate, the third value must also
        be a non-negative float if all objects in ``recrates``
        are instances of :class:`fwdpy11.Region`. However,
        if they are instead instances of :class:`fwdpy11.GeneticMapUnit`,
        then the final element in ``rates`` must be ``None``.
        See the :ref:`section <geneticmaps>` on setting
        recombination rates for details.


    .. versionadded:: 0.1.1

    .. versionchanged:: 0.2.0
        Changed this from a horrible class hierarchy
        into a much simpler, single class.

    .. versionchanged:: 0.6.0

        Updated to support :class:`fwdpy11.DiscreteDemography`

    .. versionchanged:: 0.8.0

        Refactored class internals using ``attrs``.
        Mutation and recombination rates now
        stored in :class:`fwdpy11.MutationAndRecombinationRates`
    """

    nregions = attr.ib(factory=list)
    sregions = attr.ib(factory=list)
    recregions = attr.ib(factory=list)
    rates: MutationAndRecombinationRates = attr.ib(converter=_convert_rates)
    gvalue = attr.ib(default=None)
    demography = attr.ib(
        default=fwdpy11.DiscreteDemography(), converter=_convert_demography
    )
    simlen: int = attr.ib(converter=int, default=0)
    prune_selected: bool = attr.ib(default=True)

    pself: float = attr.ib(default=0.0)  # Deprecated
    popsizes = attr.ib()  # FIXME: this is a hack from 0.6.0

    @popsizes.default
    def popsizes_default(self):
        if isinstance(self.demography, fwdpy11.DiscreteDemography):
            return None

        if isinstance(
            self.demography, fwdpy11.demographic_models.DemographicModelDetails
        ):
            return None

        # Otherwise, assume that it is a numpy array
        warnings.warn(
            "attribute popsizes is deprecated"
            " and will be removed in a future release.",
            DeprecationWarning,
        )
        return self.demography

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
        try:
            for i in value:
                attr.validators.instance_of(fwdpy11.Region)(self, attribute, i)
        except TypeError:
            try:
                for i in value:
                    attr.validators.instance_of(fwdpy11.GeneticMapUnit)(
                        self, attribute, i
                    )
            except TypeError:
                raise

    @rates.validator
    def rates_validator(self, attribute, value):
        if value.recombination_rate is None:
            for i in self.recregions:
                if not isinstance(i, fwdpy11.GeneticMapUnit):
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
            for i in value:
                attr.validators.instance_of(fwdpy11.DiploidGeneticValue)(
                    self, attribute, i
                )
        except TypeError:
            try:
                attr.validators.instance_of(fwdpy11.DiploidGeneticValue)(
                    self, attribute, value
                )
            except TypeError:
                raise

    @demography.validator
    def validate_demography(self, attribute, value):
        if isinstance(value, fwdpy11.DiscreteDemography):
            return

        if isinstance(value, fwdpy11.demographic_models.DemographicModelDetails):
            return
        raise ValueError(f"Unknown type for {attribute}: {type(value)}")

    @simlen.validator
    def validate_simlen(self, attribute, value):
        if value <= 0:
            raise ValueError(f"{attribute} must be >= 0, but we got {value} instead")

        if self.popsizes is not None and value != len(self.popsizes):
            raise ValueError(
                f"simlen must equal len(self.popsizes), but we got {value} and "
                f"{len(self.popsizes)} instead"
            )

    @pself.validator
    def validate_pself(self, attribute, value):
        if value != 0.0:
            warnings.warn(
                "attribute pself is deprecated"
                " and will be removed in a future version",
                DeprecationWarning,
            )
