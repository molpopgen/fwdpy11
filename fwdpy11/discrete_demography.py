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


@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict
@attr.s(auto_attribs=True, frozen=True, repr_ns="fwdpy11")
class MassMigration(fwdpy11._fwdpy11._ll_MassMigration):
    """
    Mass migration events.

    This class has the following attributes, whose names
    are also `kwargs` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param when: When the mass migration happens
    :type when: int
    :param source: The source deme for individuals to move
    :type source: int
    :param destination: The deme to where individuals we be moved
    :type destination: int
    :param fraction: The fraction of `source` to move to `destination`.
    :type fraction: float
    :param move_individuals: If `True`, the event moves individuals.
     If `False`, individuals are copied.
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
        warnings.warn(
            f"{self.__class__} is deprecated and will be removed soon", DeprecationWarning)
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
    warnings.warn(
        "move_individuals is deprecated and will be removed soon", DeprecationWarning)
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
    warnings.warn(
        "copy_individuals is deprecated and will be removed soon", DeprecationWarning)
    return MassMigration(when, source, destination, fraction, False, resets_growth_rate)


@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict
@attr.s(auto_attribs=True, frozen=True, repr_ns="fwdpy11")
class SetDemeSize(fwdpy11._fwdpy11._ll_SetDemeSize):
    """
    Set the size of a deme at a given time.

    This class has the following attributes, whose names
    are also `kwargs` for intitialization.  The attribute names
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
        warnings.warn(
            f"{self.__class__} is deprecated and will be removed soon", DeprecationWarning)
        super(SetDemeSize, self).__init__(
            self.when, self.deme, self.new_size, self.resets_growth_rate
        )


@attr_class_to_from_dict
@attr_add_asblack
@attr_class_pickle_with_super
@attr.s(auto_attribs=True, frozen=True, repr_ns="fwdpy11")
class SetExponentialGrowth(fwdpy11._fwdpy11._ll_SetExponentialGrowth):
    """
    Set the growth rate of a deme at a given time.

    This class has the following attributes, whose names
    are also `kwargs` for intitialization.  The attribute names
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
        warnings.warn(
            f"{self.__class__} is deprecated and will be removed soon", DeprecationWarning)
        super(SetExponentialGrowth, self).__init__(
            self.when, self.deme, self.G)


@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict
@attr.s(auto_attribs=True, frozen=True, repr_ns="fwdpy11")
class SetSelfingRate(fwdpy11._fwdpy11._ll_SetSelfingRate):
    """
    Set the selfing probability within a deme at a given time.

    This class has the following attributes, whose names
    are also `kwargs` for intitialization.  The attribute names
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
        warnings.warn(
            f"{self.__class__} is deprecated and will be removed soon", DeprecationWarning)
        super(SetSelfingRate, self).__init__(self.when, self.deme, self.S)

    def __getstate__(self):
        return self.asdict()

    def __setstate__(self, d):
        self.__dict__.update(d)
        super(SetSelfingRate, self).__init__(**d)


@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict
@attr.s(eq=False, auto_attribs=True, frozen=True, repr_ns="fwdpy11")
class MigrationMatrix(fwdpy11._fwdpy11._ll_MigrationMatrix):
    """
    The forward migration matrix for a simulation.

    This class has the following attributes, whose names
    are also `kwargs` for intitialization.  The attribute names
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
        warnings.warn(
            f"{self.__class__} is deprecated and will be removed soon", DeprecationWarning)
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
@attr.s(eq=False, auto_attribs=True, frozen=True, repr_ns="fwdpy11")
class SetMigrationRates(fwdpy11._fwdpy11._ll_SetMigrationRates):
    """
    Set the migration parameters of a simulation at a given time.
    May be used to set either the migration rates from a given
    deme or the entire migration matrix.

    This class has the following attributes, whose names
    are also `kwargs` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param when: The generation when the event occurs
    :type when: int
    :param deme: The row index of the migration matrix
    :type deme: int or None
    :param migrates: The migration rates into `deme` from all other populations.
    :type migrates: list or numpy.ndarray

    The migration rates are equivalent to the single-generation ancestry proportions
    with respect to `deme`.

    In order to change the entire migration matrix, pass `None` (or -1)
    to `deme` and a 2d `numpy.ndarray` for `migrates`.  A value of
    `None` will be converted to `-1` in this case.

    .. versionadded:: 0.5.3

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    when: int = attr.ib()
    deme: typing.Optional[int] = attr.ib(
        converter=_set_migration_rates_convert_deme)
    migrates: np.ndarray = attr.ib()

    def __attrs_post_init__(self):
        warnings.warn(
            f"{self.__class__} is deprecated and will be removed soon", DeprecationWarning)
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
            super(SetMigrationRates, self).__init__(
                self.when, self.migrates.tolist())

    def __getstate__(self):
        return self.asdict()

    def __setstate__(self, d):
        self.__dict__.update(d)
        if self.deme == -1:
            super(SetMigrationRates, self).__init__(self.when, self.migrates)
        else:
            super(SetMigrationRates, self).__init__(
                self.when, self.deme, self.migrates)

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
                    SetDemeSize(when=i, deme=0, new_size=N,
                                resets_growth_rate=True)
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
@attr.s(auto_attribs=True, frozen=True, repr_ns="fwdpy11")
class DiscreteDemography(fwdpy11._fwdpy11._ll_DiscreteDemography):
    """
    Representation of demographic events acting on
    discrete demes (sub-populations).

    This class has the following attributes, whose names
    are also `kwargs` for intitialization.  The attribute names
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
     `None`, then input are converted into an
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

    @staticmethod
    def _event_names_list() -> typing.List[str]:
        """
        Obtain the names of fields expected to contain
        lists of demographic events.
        """
        return [
            "mass_migrations",
            "set_growth_rates",
            "set_deme_sizes",
            "set_selfing_rates",
            "set_migration_rates",
        ]

    def __deepcopy__(self, memo):
        rv = self.__class__(**self.asdict())
        self._clone_state_to(rv)
        return rv

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

    def number_of_demes(self):
        demes = set()
        for e in self._timed_events():
            if e is not None:
                for i in e:
                    try:
                        demes.add(i.deme)
                    except Exception:
                        demes.add(i.source)
                        demes.add(i.destination)
        return len(demes)


def from_demes(
    dg: typing.Union[str, demes.Graph], burnin: int = 10
) -> "DemographicModelDetails":
    """
    Build a :class:`fwdpy11.DiscreteDemography` object using demes. The deme
    graph can either be a demes Graph object or a string as the filepath to a
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


@dataclass
class _DemeSize:
    deme: int
    ancestors: typing.Optional[typing.List[int]]
    when: int
    size: int
    growth_parameter: float

    def __post_init__(self):
        if self.ancestors is not None:
            self.ancestors = sorted(self.ancestors)

    def add_ancestor(self, ancestor: int):
        if self.ancestors is None:
            self.ancestors = [ancestor]
        else:
            self.ancestors.append(ancestor)
            self.ancestors = sorted(self.ancestors)


@dataclass
class _EpochData:
    deme: int
    ancestors: typing.Optional[typing.List[int]]
    initial_size: int
    growth_parameter: float

    def __post_init__(self):
        if self.ancestors is not None:
            self.ancestors = sorted(self.ancestors)

    def add_ancestor(self, ancestor: int):
        if self.ancestors is None:
            self.ancestors = [ancestor]
        else:
            self.ancestors.append(ancestor)
            self.ancestors = sorted(self.ancestors)


@dataclass
class _EventsAtTimeT:
    mass_migrations: typing.List[MassMigration]
    set_deme_sizes: typing.List[SetDemeSize]
    set_growth_rates: typing.List[SetExponentialGrowth]


@dataclass
class _SortedEvents:
    """
    Not part of the public API.

    This class defines an iterator over demographic events
    that change population sizes.
    """

    mass_migrations: typing.List[MassMigration]
    set_deme_sizes: typing.List[SetDemeSize]
    set_growth_rates: typing.List[SetExponentialGrowth]

    def __post_init__(self):
        self.mass_migrations = sorted(
            self.mass_migrations,
            # TODO: compare this to the C++ side sort function.
            key=lambda x: (x.when, x.source, x.destination,
                           not x.move_individuals),
        )

        self.set_deme_sizes = sorted(
            self.set_deme_sizes, key=lambda x: (x.when, x.deme)
        )
        self.set_growth_rates = sorted(
            self.set_growth_rates, key=lambda x: (x.when, x.deme)
        )

        for i in [self.mass_migrations, self.set_deme_sizes, self.set_growth_rates]:
            i.reverse()

        # NOTE: several properties of
        # the instance objects in these lists
        # were already checked upon initialization.

        # We do not allow things like setting the growth
        # rate > 1 time for the same deme at the same generation.
        for a, b in _pairwise(self.set_deme_sizes):
            if a.when == b.when and a.deme == b.deme:
                raise ValueError(
                    f"Multple SetDemeSize events at the same time for the same deme: {a}, {b}"
                )
        for a, b in _pairwise(self.set_growth_rates):
            if a.when == b.when and a.deme == b.deme:
                raise ValueError(
                    f"Multple SetExponentialGrowth events at the same time for the same deme: {a}, {b}"
                )
        for a, b in _pairwise(self.mass_migrations):
            # TODO: what if one is a move and the other a copy?
            if (
                a.when == b.when
                and a.source == b.source
                and a.destination == b.destination
            ):
                raise ValueError(
                    f"Multple MassMigration events at the same time for the same pair of demes: {a}, {b}"
                )

    def _get_next_event_time(self) -> typing.Optional[int]:
        times = []
        if len(self.mass_migrations) > 0:
            times.append(self.mass_migrations[-1].when)
        if len(self.set_deme_sizes) > 0:
            times.append(self.set_deme_sizes[-1].when)
        if len(self.set_growth_rates) > 0:
            times.append(self.set_growth_rates[-1].when)

        if len(times) == 0:
            return None

        return min(times)

    def _get_next_simultaneous_events(self) -> typing.Optional[_EventsAtTimeT]:
        next_event_time = self._get_next_event_time()
        if next_event_time is None:
            return None

        rv = _EventsAtTimeT([], [], [])

        for i, j in [
            (self.mass_migrations, rv.mass_migrations),
            (self.set_deme_sizes, rv.set_deme_sizes),
            (self.set_growth_rates, rv.set_growth_rates),
        ]:
            while len(i) > 0 and i[-1].when == next_event_time:
                j.append(i.pop())

        return rv

    def __iter__(self):
        return self

    def __next__(
        self,
    ):
        ne = self._get_next_simultaneous_events()
        if ne is not None:
            return ne

        raise StopIteration


def _get_N_exp_growth(N0: int, t: int, G: float) -> int:
    return int(np.rint(N0 * G**t))


def _get_epoch_end_size(
    deme: int, when: int, deme_sizes: typing.Dict[int, typing.List[_DemeSize]]
) -> int:
    if deme not in deme_sizes:
        raise ValueError(f"deme {deme} does not exist")

    assert len(deme_sizes[deme]) > 0

    return _get_N_exp_growth(
        deme_sizes[deme][-1].size,
        when - deme_sizes[deme][-1].when,
        deme_sizes[deme][-1].growth_parameter,
    )


def _pairwise(iterable: typing.Iterable) -> typing.Iterator:
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)


def _global_extinction(deme_sizes: typing.Dict[int, typing.List[_DemeSize]]):
    for _, v in deme_sizes.items():
        assert len(v) > 0
        assert v[-1].size >= 0
        if v[-1].size > 0:
            return False
    return True


# TODO / FIXME: there is so much logic duplication here.
# Need to identify patterns and abstract ops out to functions, etc..
def _process_lowlevel_mass_migrations(mass_migrations, deme_sizes):
    destination_sizes_temp = defaultdict(list)
    for mass_mig_event in mass_migrations:
        if mass_mig_event.source not in deme_sizes:
            raise ValueError(
                f"Mass migration {mass_mig_event} involves a non-existant source deme"
            )
        if deme_sizes[mass_mig_event.source][-1].size == 0:
            raise ValueError(
                f"MassMigration from an extinct source deme: {mass_mig_event}"
            )
        # This is how many individuals will mass migrate
        end_size = _get_epoch_end_size(
            mass_mig_event.source, mass_mig_event.when + 1, deme_sizes
        )
        num_source_ind = int(np.rint(end_size * mass_mig_event.fraction))
        # TODO: there's a lot of logic duplication here.
        # We may also need to check:
        # What happens if the fraction moved rounds down to 0
        # and the old growth rate was 1.0.
        # In that case, no new epoch is needed.
        if mass_mig_event.move_individuals:
            # Adjust source deme size
            if mass_mig_event.when + 1 == deme_sizes[mass_mig_event.source][-1].when:
                deme_sizes[mass_mig_event.source][-1].size -= num_source_ind
                if mass_mig_event.resets_growth_rate is True:
                    deme_sizes[mass_mig_event.source][-1].growth_parameter = 1.0
            else:
                new_size = end_size - num_source_ind
                if new_size < 0:
                    raise ValueError(
                        f"MassMigration results in a source deme size < 0, {mass_mig_event}"
                    )
                if mass_mig_event.resets_growth_rate is True:
                    G = 1.0
                else:
                    G = deme_sizes[mass_mig_event.source][-1].growth_parameter

                if deme_sizes[mass_mig_event.source][-1].ancestors is not None:
                    new_ancestors = deme_sizes[mass_mig_event.source][
                        -1
                    ].ancestors.copy()
                    if mass_mig_event.source not in new_ancestors:
                        new_ancestors.append(mass_mig_event.source)
                        new_ancestors = sorted(new_ancestors)
                else:
                    new_ancestors = [mass_mig_event.source]

                if mass_mig_event.fraction == 1.0:
                    new_ancestors = None

                deme_sizes[mass_mig_event.source].append(
                    _DemeSize(
                        mass_mig_event.source,
                        new_ancestors,
                        mass_mig_event.when + 1,
                        new_size,
                        G,
                    )
                )
        else:  # A copy
            if mass_mig_event.when + 1 == deme_sizes[mass_mig_event.source][-1].when:
                if mass_mig_event.resets_growth_rate is True:
                    deme_sizes[mass_mig_event.source][-1].growth_parameter = 1.0
            elif deme_sizes[mass_mig_event.source][-1].growth_parameter != 1.0:
                # Only add a new epoch if the growth rate actually
                # would be reset
                new_size = end_size
                if mass_mig_event.resets_growth_rate is True:
                    G = 1.0
                else:
                    G = deme_sizes[mass_mig_event.source][-1].growth_parameter

                if deme_sizes[mass_mig_event.source][-1].ancestors is not None:
                    new_ancestors = deme_sizes[mass_mig_event.source][
                        -1
                    ].ancestors.copy()
                    if mass_mig_event.source not in new_ancestors:
                        new_ancestors.append(mass_mig_event.source)
                        new_ancestors = sorted(new_ancestors)
                else:
                    new_ancestors = [mass_mig_event.source]
                deme_sizes[mass_mig_event.source].append(
                    _DemeSize(
                        mass_mig_event.source,
                        new_ancestors,
                        mass_mig_event.when + 1,
                        new_size,
                        G,
                    )
                )

        if mass_mig_event.destination not in deme_sizes:
            # This is a new deme, so we initialize
            # it with growth parameter of 1.0
            destination_sizes_temp[mass_mig_event.destination].append(
                _DemeSize(
                    mass_mig_event.destination,
                    [mass_mig_event.source],  # A new deme
                    mass_mig_event.when + 1,
                    num_source_ind,
                    1.0,
                )
            )
        else:
            assert len(deme_sizes[mass_mig_event.destination]) > 0

            if mass_mig_event.destination in destination_sizes_temp:
                # There are other mass migrations into this deme at this time.
                assert (
                    mass_mig_event.when + 1
                    == destination_sizes_temp[mass_mig_event.destination][-1].when
                )
                assert (
                    destination_sizes_temp[mass_mig_event.destination][-1].ancestors
                    is not None
                )
                end_size = _get_epoch_end_size(
                    mass_mig_event.destination,
                    mass_mig_event.when + 1,
                    destination_sizes_temp,
                )
                new_size = (
                    destination_sizes_temp[mass_mig_event.destination][-1].size
                    + num_source_ind
                )
                if new_size < 0:
                    raise ValueError(
                        f"MassMigration results in a source deme size < 0, {mass_mig_event}, {new_size}"
                    )
                if mass_mig_event.resets_growth_rate is True:
                    G = 1.0
                else:
                    G = destination_sizes_temp[mass_mig_event.destination][
                        -1
                    ].growth_parameter
                if (
                    destination_sizes_temp[mass_mig_event.destination][-1].ancestors
                    is not None
                ):
                    new_ancestors = destination_sizes_temp[mass_mig_event.destination][
                        -1
                    ].ancestors.copy()
                    if mass_mig_event.source not in new_ancestors:
                        new_ancestors.append(mass_mig_event.source)
                        new_ancestors = sorted(new_ancestors)
                else:
                    new_ancestors = [mass_mig_event.source]

                last = destination_sizes_temp[mass_mig_event.destination][-1]
                last.size = new_size
                last.growth_parameter = G
                last.ancestors = new_ancestors
            else:
                # This is the first mass migration at this time
                # into this deme
                end_size = _get_epoch_end_size(
                    mass_mig_event.destination,
                    mass_mig_event.when + 1,
                    deme_sizes,
                )
                new_size = end_size + num_source_ind
                if new_size < 0:
                    raise ValueError(
                        f"MassMigration results in a source deme size < 0, {mass_mig_event}, {new_size}"
                    )
                if mass_mig_event.resets_growth_rate is True:
                    G = 1.0
                else:
                    G = deme_sizes[mass_mig_event.destination][-1].growth_parameter
                if deme_sizes[mass_mig_event.destination][-1].ancestors is not None:
                    new_ancestors = sorted(
                        [mass_mig_event.destination, mass_mig_event.source]
                    )
                else:
                    new_ancestors = [mass_mig_event.source]
                destination_sizes_temp[mass_mig_event.destination].append(
                    _DemeSize(
                        mass_mig_event.destination,
                        new_ancestors,
                        mass_mig_event.when + 1,
                        new_size,
                        G,
                    )
                )

    for k, v in destination_sizes_temp.items():
        deme_sizes[k].extend(v)

    return len(mass_migrations)


def _process_lowlevel_set_deme_sizes(set_deme_sizes, deme_sizes):
    for event in set_deme_sizes:
        if event.resets_growth_rate is True and event.deme in deme_sizes:
            assert len(deme_sizes[event.deme]) > 0
            if event.when + 1 != deme_sizes[event.deme][-1].when:
                if deme_sizes[event.deme][-1].ancestors is None:
                    new_ancestors = None
                elif deme_sizes[event.deme][-1].size == 0:
                    new_ancestors = None
                else:
                    new_ancestors = [event.deme]

                deme_sizes[event.deme].append(
                    _DemeSize(
                        event.deme,
                        new_ancestors,
                        event.when + 1,
                        event.new_size,
                        1.0,
                    )
                )
            else:
                deme_sizes[event.deme][-1].size = event.new_size
                deme_sizes[event.deme][-1].growth_parameter = 1.0
                if deme_sizes[event.deme][-1].size == 0:
                    deme_sizes[event.deme][-1].ancestors = None
        else:
            # Inherit the growth rate from the previous epoch
            if event.deme in deme_sizes:
                assert len(deme_sizes[event.deme]) > 0
                lastG = deme_sizes[event.deme][-1].growth_parameter
                if event.when + 1 != deme_sizes[event.deme][-1].when:
                    if deme_sizes[event.deme][-1].ancestors is None:
                        new_ancestors = None
                    elif deme_sizes[event.deme][-1].size == 0:
                        new_ancestors = None
                    else:
                        new_ancestors = [event.deme]

                    deme_sizes[event.deme].append(
                        _DemeSize(
                            event.deme,
                            new_ancestors,
                            event.when + 1,
                            event.new_size,
                            lastG,
                        )
                    )
                else:
                    deme_sizes[event.deme][-1].size = event.new_size
                    deme_sizes[event.deme][-1].growth_parameter = lastG
                    if event.new_size == 0:
                        deme_sizes[event.deme][-1].ancestors = None
            else:
                # There is no previous epoch for this deme,
                # so we initialize with G = 1.0
                deme_sizes[event.deme].append(
                    _DemeSize(event.deme, None, event.when +
                              1, event.new_size, 1.0)
                )

    return len(set_deme_sizes)


def _process_lowlevel_set_growth_rates(set_growth_rates, deme_sizes):
    for event in set_growth_rates:
        if (
            event.deme not in deme_sizes
            or len(deme_sizes[event.deme]) == 0
            or deme_sizes[event.deme][-1].size == 0
        ):
            raise ValueError(
                f"Event {event} being applied to a non-existant deme")
        if event.when + 1 == deme_sizes[event.deme][-1].when:
            # the rate is concurrent w/another event, so
            # just change the value
            deme_sizes[event.deme][-1].growth_parameter = event.G
        else:
            assert event.when + 1 > deme_sizes[event.deme][-1].when
            N0 = _get_N_exp_growth(
                deme_sizes[event.deme][-1].size,
                event.when + 1 - deme_sizes[event.deme][-1].when,
                deme_sizes[event.deme][-1].growth_parameter,
            )
            deme_sizes[event.deme].append(
                _DemeSize(event.deme, [event.deme],
                          event.when + 1, N0, event.G)
            )

    return len(set_growth_rates)


@dataclass
class _DemeSizeHistory:
    """
    This class is not part of the public API.

    It is used by various parts of the Python back end
    for validating demographic models.
    """

    epochs: intervaltree.IntervalTree
    model_times: typing.Optional[typing.ForwardRef("_ModelTimes")]

    def __post_init__(self):
        self._validate_epoch_ancestry()

    @staticmethod
    def from_demes_graph(
        g: demes.Graph,
        burnin_generation: int,
        idmap: typing.Dict,
        model_times: typing.ForwardRef("_ModelTimes"),
    ):
        """
        NOTE: this function provides a limited data structure
              that only allows existence querying at given times.

        """
        itree = intervaltree.IntervalTree()
        for d in g.demes:
            for e in d.epochs:
                start = model_times.convert_time(e.start_time)
                # NOTE: this is to ensure that the interval tree
                # correctly contains the 1/2 open interval
                # during which individuals exist in this deme
                if start > 0:
                    start = start + 1
                end = model_times.convert_time(e.end_time) + 1
                if not end > start:
                    raise ValueError(
                        f"epoch end time must be > start, got end: {end}, start: {start}"
                    )
                # NOTE: 1.0 is a HACK for growth rate and is WRONG
                itree[start:end] = _EpochData(
                    idmap[d.name], None, int(e.start_size), 1.0
                )

        rv = _DemeSizeHistory(itree, model_times)
        return rv

        # raise NotImplementedError(
        #     "_DemeSizeHistory.from_demes_graph is not implemented"
        # )

    # TODO:
    # * We need to track the ancestor deme of each epoch.
    #   Without this, we cannot distinguish a mass migration
    #   from a deme size change w/no accompanying migration
    #   to provide ancestry.
    # * We've done some of the above, but probably incorrectly:
    #   We need tests of 3- and 4- way mass-migration.
    #   More generally, we need to account for multi-way ancestry
    @staticmethod
    def from_lowlevel(
        initial_sizes: typing.Dict[int, int],
        *,
        mass_migrations: typing.Optional[typing.List[MassMigration]] = None,
        set_deme_sizes: typing.Optional[typing.List[SetDemeSize]] = None,
        set_growth_rates: typing.Optional[typing.List[SetExponentialGrowth]] = None,
        total_simulation_length: typing.Optional[int] = None,
    ):
        # TODO: we repeat ourselves a lot
        # in how we access deme_sizes.
        # This should probably be abstracted out
        # to a class containing a dict.
        deme_sizes = defaultdict(list)
        for deme, size in initial_sizes.items():
            if deme < 0:
                raise ValueError(
                    f"all deme indexes must be non-negative, got {deme}")
            if size <= 0:
                raise ValueError(f"deme {deme} size must be > 0, got {size}")
            deme_sizes[deme].append(_DemeSize(deme, [deme], 0, size, 1.0))

        events = _SortedEvents(
            mass_migrations if mass_migrations is not None else [],
            set_deme_sizes if set_deme_sizes is not None else [],
            set_growth_rates if set_growth_rates is not None else [],
        )

        num_events_processed = 0
        global_extinction = False
        num_events_expected = (
            len(events.mass_migrations)
            + len(events.set_growth_rates)
            + len(events.set_deme_sizes)
        )
        for next_events in events:
            global_extinction = _global_extinction(deme_sizes)
            if global_extinction is True:
                break

            num_events_processed += _process_lowlevel_mass_migrations(
                next_events.mass_migrations, deme_sizes
            )
            num_events_processed += _process_lowlevel_set_deme_sizes(
                next_events.set_deme_sizes, deme_sizes
            )
            num_events_processed += _process_lowlevel_set_growth_rates(
                next_events.set_growth_rates, deme_sizes
            )

        if not global_extinction:
            assert (
                num_events_processed == num_events_expected
            ), f"{num_events_processed} {num_events_expected} {events}"

        itree = intervaltree.IntervalTree()

        for k, v in deme_sizes.items():
            assert len(v) > 0
            for a, b in _pairwise(v):
                itree[a.when: b.when] = _EpochData(
                    k, a.ancestors, a.size, a.growth_parameter
                )
            if total_simulation_length is not None:
                if v[-1].when < total_simulation_length and v[-1].size > 0:
                    itree[v[-1].when: total_simulation_length + 1] = _EpochData(
                        k, v[-1].ancestors, v[-1].size, v[-1].growth_parameter
                    )
            else:
                if v[-1].size > 0:
                    itree[v[-1].when: int(np.iinfo(np.uint32).max)] = _EpochData(
                        k, v[-1].ancestors, v[-1].size, v[-1].growth_parameter
                    )
        rv = _DemeSizeHistory(itree, None)
        return rv

    def _validate_epoch_ancestry(self):
        """
        If this function fails, the most likely
        cause is a logic error in setting up the
        ancestors.
        """
        for interval in self.epochs:
            if interval.data.ancestors is not None:
                for ancestor in interval.data.ancestors:
                    if interval.begin == 0:
                        if interval.data.deme != ancestor:
                            raise RuntimeError("epoch ancestry error")
                        w = interval.begin
                    else:
                        w = interval.begin - 1
                    if not self.deme_exists_at(ancestor, w):
                        raise RuntimeError("epoch ancestry error")

    def demes_at(self, generation: int):
        if generation < 0:
            raise ValueError("generation must be non-negative")
        for interval in self.epochs[generation]:
            yield interval.data.deme

    def deme_exists_at(self, deme: int, generation: int) -> bool:
        if deme < 0:
            raise ValueError("deme must be non-negative")
        for d in self.demes_at(generation):
            if d == deme:
                return True
        return False

    def deme_size_at(self, deme: int, generation: int) -> typing.Optional[int]:
        if not self.deme_exists_at(deme, generation):
            return None

        for interval in self.epochs[generation]:
            if interval.data.deme == deme:
                data = interval.data
                t = generation - interval.begin + 1
                return int(np.rint(data.initial_size * data.growth_parameter**t))

    def demes_coexist_at(self, demes: typing.Tuple[int, int], generation: int) -> bool:
        for deme in demes:
            if deme < 0:
                raise ValueError("deme must be non-negative")
        temp = [i for i in self.demes_at(generation)]
        return all([i in temp for i in demes])
