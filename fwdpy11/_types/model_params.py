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


@attr.s(kw_only=True, frozen=True)
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
    :type rates: list
    :param demography: The demographic model to simulate
    :type demography: fwdpy11.DiscreteDemography
    :param simlen: The number of time steps to evolve
    :type simlen: int
    :param prune_selected: If ``True``, remove selected fixations from
                           the population when they are first detected
                           as fixed.
    :type prune_selected: bool

    The following attributes exist, and may be set as ``kwargs``,
    but are best set via ``rates``.  Currently, these attributes
    are only used internally.  As such, they are prone to
    refactoring in future versions.

    :param mutrate_n: The neutral mutation rate
    :type mutrate_n: float
    :param mutrate_s: The neutral mutation rate
    :type mutrate_s: float
    :param recrate: The recombination rate
    :type recrate: float or None

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

        The ``rates`` field must be a list of length three (3).
        The first two values must be non-negative floats. For
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
    """

    nregions = attr.ib(factory=list)
    sregions = attr.ib(factory=list)
    recregions = attr.ib(factory=list)
    rates = attr.ib(factory=list, converter=list)
    gvalue = attr.ib(default=None)
    demography = attr.ib(default=fwdpy11.DiscreteDemography())
    simlen: int = attr.ib(converter=int)
    prune_selected: bool = attr.ib(default=True)

    pself: float = attr.ib(default=0.0)  # Deprecated
    popsizes = attr.ib()  # FIXME: this is a hack from 0.6.0

    # These attributes are the elements
    # store in rates
    mutrate_n: float = attr.ib(
        converter=float, validator=attr.validators.instance_of(float)
    )
    mutrate_s: float = attr.ib(
        converter=float, validator=attr.validators.instance_of(float)
    )
    recrate = attr.ib()

    @mutrate_n.default
    def mutrate_n_default(self):
        return self.rates[0]

    @mutrate_s.default
    def mutrate_s_default(self):
        return self.rates[1]

    @recrate.default
    def recrate_default(self):
        if self.rates[2] is None:
            return None
        return float(self.rates[2])

    @popsizes.default
    def popsizes_default(self):
        if isinstance(self.demography, fwdpy11.DiscreteDemography):
            return None

        # Otherwise, assume that it is a numpy array
        warnings.warn(
            "attribute popsizes is being considered for"
            " deprecation (along with simulations without tree sequences)",
            PendingDeprecationWarning,
        )
        return self.demography

    @simlen.default
    def simlen_default(self):
        if isinstance(self.demography, fwdpy11.DiscreteDemography):
            return 0
        return len(self.demography)

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
        if len(value) != 3:
            raise ValueError(f"rates must be of length 3, but we got {value} instead")

    @mutrate_n.validator
    @mutrate_s.validator
    def individual_rate_validator(self, attribute, value):
        # NOTE: not decorated as recrate validator b/c that
        # attribute is allowed to be None
        if not np.isfinite(value):
            raise ValueError(f"{attribute} must be finite, but we got {value} instead")
        if value < 0.0:
            raise ValueError(f"{attribute} must be >= 0.0, but we got {value} instead")

    @mutrate_n.validator
    def validate_mutrate_n(self, attribute, value):
        if value > 0 and len(self.nregions) == 0:
            raise ValueError(
                "mutation rate to neutral variants is > 0 but no Regions are defined"
            )

    @mutrate_s.validator
    def validate_mutrate_s(self, attribute, value):
        if value > 0 and len(self.sregions) == 0:
            raise ValueError(
                "mutation rate to selected variants is > 0 but no Sregions are defined"
            )

    @recrate.validator
    def validate_recrate(self, attribute, value):
        if value is None:
            for i in self.recregions:
                if not isinstance(i, fwdpy11.GeneticMapUnit):
                    raise ValueError(
                        f"recrate of {value} must be paired with"
                        " instances of fwdpy11.GeneticMapUnit"
                    )
        else:
            for i in self.recregions:
                if not isinstance(i, fwdpy11.Region):
                    raise ValueError(
                        f"recrate of {value} must be paired with"
                        " instances of fwdpy11.Region"
                    )
            attr.validators.instance_of(float)(self, attribute, value)
            self.individual_rate_validator(attribute, value)

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
                "attribute pself is being considered for"
                " deprecation (along with simulations without tree sequences)",
                PendingDeprecationWarning,
            )

    def as_dict(self):
        """
        Return instance attributes as a :class:`dict`.

        .. versionadded:: 0.8.0
        """
        return attr.asdict(self)
