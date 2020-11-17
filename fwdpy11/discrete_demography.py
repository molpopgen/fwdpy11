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

import typing

import attr
import numpy as np

import fwdpy11

from .class_decorators import (attr_add_asblack, attr_class_pickle_with_super,
                               attr_class_to_from_dict,
                               attr_class_to_from_dict_no_recurse)

_common_attr_attribs = {"frozen": True, "auto_attribs": True, "repr_ns": "fwdpy11"}


@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict
@attr.s(**_common_attr_attribs)
class MassMigration(fwdpy11._fwdpy11._ll_MassMigration):
    """
    Mass migration events.

    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param when: When the mass migration happens
    :type when: int
    :param source: The source deme for individuals to move
    :type source: int
    :param destination: The deme to where individuals we be moved
    :type destination: int
    :param fraction: The fraction of `source` to move to `destination`.
    :type fraction: float
    :param move_individuals: If ``True``, the event moves individuals.
                             If ``False``, individuals are copied.
    :type move_individuals: bool
    :param resets_growth_rate: (True) Whether or not to reset the
                               growth rates of `source` and `destination` to
                               :data:`fwdpy11.NOGROWTH`
    :type resets_growth_rate: bool

    .. note::

        It is a bit simpler to construct instances of this class
        using :func:`fwdpy11.move_individuals` or :func:`fwdpy11.copy_individuals`.

    .. versionadded:: 0.5.3

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    when: int = attr.ib()
    source: int = attr.ib()
    destination: int = attr.ib()
    fraction: float = attr.ib()
    move_individuals: bool = attr.ib(default=True)
    resets_growth_rate: bool = attr.ib(default=True)

    def __attrs_post_init__(self):
        super(MassMigration, self).__init__(
            self.when,
            self.source,
            self.destination,
            self.fraction,
            self.move_individuals,
            self.resets_growth_rate,
        )


def move_individuals(
    when: int,
    source: int,
    destination: int,
    fraction: float,
    resets_growth_rate: bool = True,
) -> MassMigration:
    """
    :param when: The generation when the event occurs
    :type when: int
    :param source: The source deme for individuals to move
    :type source: int
    :param destination: The deme to where individuals we be moved
    :type destination: int
    :param fraction: The fraction of `source` to move to `destination`.
    :type fraction: float
    :param resets_growth_rate: (True) Whether or not to reset the
                               growth rates of `source` and `destination` to
                               :data:`fwdpy11.NOGROWTH`
    :type resets_growth_rate: bool

    :rtype: :class:`fwdpy11.MassMigration`
    """
    return MassMigration(when, source, destination, fraction, True, resets_growth_rate)


def copy_individuals(
    when: int,
    source: int,
    destination: int,
    fraction: float,
    resets_growth_rate: bool = True,
) -> MassMigration:
    """
    :param when: The generation when the event occurs
    :type when: int
    :param source: The source deme for individuals to copy
    :type source: int
    :param destination: The deme to where individuals we be copied
    :type destination: int
    :param fraction: The fraction of `source` to copy to `destination`.
    :type fraction: float
    :param resets_growth_rate: (True) Whether or not to reset the
                               growth rates of `source` and `destination` to
                               :data:`fwdpy11.NOGROWTH`
    :type resets_growth_rate: bool

    :rtype: :class:`fwdpy11.MassMigration`
    """
    return MassMigration(when, source, destination, fraction, False, resets_growth_rate)


@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict
@attr.s(**_common_attr_attribs)
class SetDemeSize(fwdpy11._fwdpy11._ll_SetDemeSize):
    """
    Set the size of a deme at a given time.

    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param when: The generation when the event occurs
    :type when: int
    :param deme: The deme whose size will change
    :type deme: int
    :param new_size: The new size
    :type new_size: int
    :param resets_growth_rate: (True) If deme size change resets
                               growth rate to :data:`fwdpy11.NOGROWTH`
    :type resets_growth_rate: bool

    .. versionadded:: 0.5.3

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    when: int = attr.ib()
    deme: int = attr.ib()
    new_size: int = attr.ib()
    resets_growth_rate: bool = attr.ib(default=True)

    def __attrs_post_init__(self):
        super(SetDemeSize, self).__init__(
            self.when, self.deme, self.new_size, self.resets_growth_rate
        )


@attr_class_to_from_dict
@attr_add_asblack
@attr_class_pickle_with_super
@attr.s(**_common_attr_attribs)
class SetExponentialGrowth(fwdpy11._fwdpy11._ll_SetExponentialGrowth):
    """
    Set the growth rate of a deme at a given time.

    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param when: The generation when the event occurs
    :type when: int
    :param deme: The deme whose growth rate will change
    :type deme: int
    :param G: The new growth rate
    :type G: float

    .. versionadded:: 0.5.3

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    when: int = attr.ib()
    deme: int = attr.ib()
    G: float = attr.ib()

    def __attrs_post_init__(self):
        super(SetExponentialGrowth, self).__init__(self.when, self.deme, self.G)


@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict
@attr.s(**_common_attr_attribs)
class SetSelfingRate(fwdpy11._fwdpy11._ll_SetSelfingRate):
    """
    Set the selfing probability within a deme at a given time.

    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param when: The generation when the event occurs
    :type when: int
    :param deme: The deme whose selfing probability will change
    :type deme: int
    :param S: The new selfing probability
    :type S: float

    .. versionadded:: 0.5.3

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    when: int = attr.ib()
    deme: int = attr.ib()
    S: float = attr.ib()

    def __attrs_post_init__(self):
        super(SetSelfingRate, self).__init__(self.when, self.deme, self.S)

    def __getstate__(self):
        return self.asdict()

    def __setstate__(self, d):
        self.__dict__.update(d)
        super(SetSelfingRate, self).__init__(**d)


@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict
@attr.s(eq=False, **_common_attr_attribs)
class MigrationMatrix(fwdpy11._fwdpy11._ll_MigrationMatrix):
    """
    The forward migration matrix for a simulation.

    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param migmatrix: A square matrix of non-negative floats.
    :type migmatrix: numpy.ndarray
    :param scaled: (True) If entries in `migmatrix` will be
                   multiplied by deme sizes during simulation
    :type scaled: bool

    .. versionadded:: 0.5.3

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    migmatrix: np.ndarray = attr.ib()
    scaled: bool = attr.ib(default=False)

    def __attrs_post_init__(self):
        super(MigrationMatrix, self).__init__(self.migmatrix, self.scaled)

    def __eq__(self, other):
        return self.scaled == other.scaled and np.array_equal(
            self.migmatrix, other.migmatrix
        )

    @property
    def M(self):
        return self.migmatrix

    @property
    def shape(self):
        return self.migmatrix.shape


def _set_migration_rates_convert_deme(i: typing.Optional[int]) -> int:
    if i is None:
        return -1
    return int(i)


@attr_add_asblack
@attr_class_to_from_dict
@attr.s(eq=False, **_common_attr_attribs)
class SetMigrationRates(fwdpy11._fwdpy11._ll_SetMigrationRates):
    """
    Set the migration parameters of a simulation at a given time.
    May be used to set either the migration rates from a given
    deme or the entire migration matrix.

    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param when: The generation when the event occurs
    :type when: int
    :param deme: The row index of the migration matrix
    :type deme: int or None
    :param migrates: The migration rates into `deme` from all other populations.
    :type migrates: list or numpy.ndarray

    The migration rates are equivalent to the single-generation ancestry proportions
    with respect to ``deme``.

    In order to change the entire migration matrix, pass ``None`` (or -1)
    to ``deme`` and a 2d ``numpy.ndarray`` for ``migrates``.  A value of
    ``None`` will be converted to ``-1`` in this case.

    .. versionadded:: 0.5.3

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    when: int = attr.ib()
    deme: typing.Optional[int] = attr.ib(converter=_set_migration_rates_convert_deme)
    migrates: np.ndarray = attr.ib()

    def __attrs_post_init__(self):
        if self.deme >= 0:
            try:
                super(SetMigrationRates, self).__init__(
                    self.when, self.deme, self.migrates.tolist()
                )
            except AttributeError:
                super(SetMigrationRates, self).__init__(
                    self.when, self.deme, self.migrates
                )
        else:
            super(SetMigrationRates, self).__init__(self.when, self.migrates.tolist())

    def __getstate__(self):
        return self.asdict()

    def __setstate__(self, d):
        self.__dict__.update(d)
        if self.deme == -1:
            super(SetMigrationRates, self).__init__(self.when, self.migrates)
        else:
            super(SetMigrationRates, self).__init__(self.when, self.deme, self.migrates)

    def __eq__(self, other):
        return (
            self.deme == other.deme
            and self.when == other.when
            and np.array_equal(self.migrates, other.migrates)
        )


def _convert_demographic_events_list(
    o: typing.Optional[typing.List],
) -> typing.Optional[typing.List]:
    if o is None:
        return o

    try:
        rv = sorted(o, key=lambda x: (x.when, x.move_individuals))
    except AttributeError:
        rv = sorted(o, key=lambda x: (x.when))

    return rv


def _convert_set_deme_sizes(
    o: typing.Optional[np.ndarray],
) -> typing.Optional[np.ndarray]:
    if o is None:
        return o
    if isinstance(o, np.ndarray):
        N = o[0]
        rv = [SetDemeSize(when=0, deme=0, new_size=N, resets_growth_rate=True)]
        for i in range(1, len(o)):
            if o[i] != N:
                N = o[i]
                rv.append(
                    SetDemeSize(when=i, deme=0, new_size=N, resets_growth_rate=True)
                )

        return rv

    return _convert_demographic_events_list(o)


def _convert_migmatrix(o):
    if o is None:
        return o
    try:
        return MigrationMatrix(o.migmatrix, o.scaled)
    except AttributeError:
        pass
    if isinstance(o, tuple):
        return MigrationMatrix(*o)
    return MigrationMatrix(o)


@attr_add_asblack
@attr_class_to_from_dict_no_recurse
@attr.s(**_common_attr_attribs)
class DiscreteDemography(fwdpy11._fwdpy11._ll_DiscreteDemography):
    """
    Representation of demographic events acting on
    discrete demes (sub-populations).

    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param mass_migrations: Instances of :class:`fwdpy11.MassMigration`
    :type mass_migrations: None or list[fwdpy11.MassMigration]
    :param set_growth_rates: Instances of :class:`fwdpy11.SetExponentialGrowth`
    :type set_growth_rates: None or list[fwdpy11.SetExponentialGrowth]
    :param set_deme_sizes: Instances of :class:`fwdpy11.SetDemeSize`
    :type set_deme_sizes: None or list[fwdpy11.SetDemeSize]
    :param set_selfing_rates: Instances of :class:`fwdpy11.SetSelfingRate`
    :type set_selfing_rates: None or list[fwdpy11.SetSelfingRate]
    :param migmatrix: A migraton matrix. See :ref:`migration`. If not
                      ``None``, then input are converted into an
                      instance of :class:`fwdpy11.MigrationMatrix`.
    :type migmatrix: None or numpy.ndarray or fwdpy11.MigrationMatrix
    :param set_migration_rates: Instances of :class:`fwdpy11.SetMigrationRates`
    :type set_migration_rates: None or list[fwdpy11.SetMigrationRates]

    .. versionadded:: 0.5.3

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    mass_migrations: typing.Optional[typing.List[MassMigration]] = attr.ib(
        default=None, converter=_convert_demographic_events_list
    )
    set_growth_rates: typing.Optional[typing.List[SetExponentialGrowth]] = attr.ib(
        default=None, converter=_convert_demographic_events_list
    )
    set_deme_sizes: typing.Optional[typing.List[SetDemeSize]] = attr.ib(
        default=None, converter=_convert_set_deme_sizes
    )
    set_selfing_rates: typing.Optional[typing.List[SetSelfingRate]] = attr.ib(
        default=None, converter=_convert_demographic_events_list
    )
    migmatrix: typing.Optional[typing.Union[np.ndarray, MigrationMatrix]] = attr.ib(
        default=None, converter=_convert_migmatrix
    )
    set_migration_rates: typing.Optional[typing.List[SetMigrationRates]] = attr.ib(
        default=None, converter=_convert_demographic_events_list
    )

    def __attrs_post_init__(self):
        super(DiscreteDemography, self).__init__(
            mass_migrations=self.mass_migrations,
            set_growth_rates=self.set_growth_rates,
            set_deme_sizes=self.set_deme_sizes,
            set_selfing_rates=self.set_selfing_rates,
            migmatrix=self.migmatrix,
            set_migration_rates=self.set_migration_rates,
        )

    def __getstate__(self):
        return (self.asdict(), self._state_asdict())

    def __setstate__(self, t):
        self.__dict__.update(t[0])
        super(DiscreteDemography, self).__init__(**t[0])
        self._reset_state(t[1])

    def _timed_events(self):
        for i in [
            self.mass_migrations,
            self.set_deme_sizes,
            self.set_growth_rates,
            self.set_selfing_rates,
            self.set_migration_rates,
        ]:
            yield i
