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
import collections
import copy
import typing
import warnings

import attr
import numpy as np

import fwdpy11


def _create_event_list(o):
    try:
        d = o.model.asdict()
        dc = copy.deepcopy(d)
        return fwdpy11.DiscreteDemography(**dc)
    except AttributeError:
        d = o.asdict()
        dc = copy.deepcopy(d)
        return fwdpy11.DiscreteDemography(**dc)


def _create_initial_deme_sizes(o):
    try:
        md = np.array(o.diploid_metadata, copy=False)
        return np.unique(md["deme"], return_counts=True)[1].tolist()
    except AttributeError:
        return o


def _convert_simlen(simlen):
    if simlen is None:
        return simlen
    return int(simlen)


@attr.s()
class DemographyDebugger(object):
    """
    Efficiently debug demographic events.

    The class executes a mock simulation
    of demographic events in order to quickly
    detect errors. It also generates a "report"
    of the model details for double-checking.

    If model errors are found, a ``ValueError`` exception
    is raised. The goal is to catch the types of errors that would
    result in an exception during a simulation.

    The class is initialized with the following arguments,
    which may be either positional or keyword:

    :param initial_deme_sizes: The initial sizes of each deme
    :type initial_deme_sizes: list or fwdpy11.DiploidPopulation
    :param events: The demographic model
    :type events: fwdpy11.DiscreteDemography or
                  fwdpy11.demographic_models.DemographicModelDetails
    :param simlen: The length of the simulation.  Defaults to `None`.
    :type simlen: int
    :param deme_labels: A map from deme index to a printable name
    :type deme_labels: dict

    .. versionadded:: 0.6.0

    .. versionchanged:: 0.8.1

        Initialization now done with attr.
        A list of initial deme sizes is now accepted.
    """

    initial_deme_sizes = attr.ib(converter=_create_initial_deme_sizes)
    _events = attr.ib(converter=_create_event_list)
    simlen = attr.ib(converter=_convert_simlen, default=None)
    deme_labels = attr.ib(default=None)

    def __attrs_post_init__(self):
        self._validate_migration_rate_change_lengths(self._events)

        self.maxdemes = self._get_maxdemes()

        if self.maxdemes < 1:
            raise ValueError(
                "Invalid number of " "demes in simulation: {}".format(self.maxdemes)
            )

        self.current_deme_sizes = np.zeros(self.maxdemes, dtype=np.uint32)
        for i, j in enumerate(self.initial_deme_sizes):
            self.current_deme_sizes[i] = j

        self.growth_rates = np.array([fwdpy11.NOGROWTH] * self.maxdemes)
        self.growth_onset_times = np.zeros(self.maxdemes, dtype=np.uint32)
        self.selfing_rates = np.zeros(self.maxdemes)
        self.growth_initial_sizes = np.copy(self.current_deme_sizes)
        # NOTE/FIXME: some of the logic re: went_extinct and has_metadata
        # is redundant, convoluted, or both.
        self.went_extinct = np.zeros(self.maxdemes, dtype=np.int32)
        self.has_metadata = np.zeros(self.maxdemes, dtype=np.int32)
        self.has_metadata[(self.current_deme_sizes > 0)] = 1

        # The real work
        self._report = None
        self._process_demographic_model(self._events, self.simlen)

    @initial_deme_sizes.validator
    def validate_initial_deme_sizes(self, attribute, value):
        if len(value) == 0:
            raise ValueError("Initial deme sizes not provided")
        if any([i < 0 for i in value]) is True:
            raise ValueError("Initial deme sizes cannot be negative")
        for i in value:
            attr.validators.instance_of(int)(self, attribute, i)

    @simlen.validator
    def validate_simlen(self, attribute, value):
        if value is not None:
            attr.validators.instance_of(int)(self, attribute, value)
            if value <= 0:
                raise ValueError("simlen must be None or > 0")

    def _get_maxdemes(self):
        """
        The maximum number of demes the sim can ever see.
        """
        max_from_md = len(self.initial_deme_sizes)
        max_from_events = -1

        def update_max_from_events(m, e):
            if e is None:
                return m
            for i in e:
                try:
                    m = max(m, i.deme)
                except AttributeError:
                    m = max(m, i.source)
                    m = max(m, i.destination)
            return m

        for i in [
            self._events.mass_migrations,
            self._events.set_growth_rates,
            self._events.set_deme_sizes,
            self._events.set_selfing_rates,
            self._events.set_migration_rates,
        ]:
            max_from_events = update_max_from_events(max_from_events, i)

        current_max = max(max_from_md, max_from_events) + 1
        if self._events.migmatrix is None:
            return current_max

        if self._events.migmatrix.shape[0] < current_max:
            raise ValueError(
                f"The MigrationMatrix shape, {self._events.migmatrix.shape}, "
                f"does not match the max number of "
                f"demes present in the "
                f"simulation, {current_max}"
            )

        return max(current_max, self._events.migmatrix.shape[0])

    def _validate_migration_rate_change_lengths(self, events):
        """
        Various checks on the lengths of new migration rates.
        """
        if (
            self._events.migmatrix is None
            and events.set_migration_rates is not None
            and len(events.set_migration_rates) > 0
        ):
            raise ValueError(
                "migration rate changes are " "set but there is no migration matrix"
            )
        if events.set_migration_rates is not None:
            for i in events.set_migration_rates:
                if i.deme >= 0:
                    if len(i.migrates) != self._events.migmatrix.shape[0]:
                        raise ValueError("Migration rates mismatch")
                else:  # Are replacing the entire matrix
                    if len(i.migrates.flatten()) != len(
                        self._events.migmatrix.M.flatten()
                    ):
                        raise ValueError("Migration rates mismatch")

    def _get_event_names(self, events):
        # NOTE: we rely on the fact that all public fields of
        # DiscreteDemography are event lists or the migration
        # matrix.  Thus, we filter out migmatrix, all
        # magic fxns, and some others that aren't relevant.
        # FIXME: a more Pythonic approach would be only to
        # pull out attributes that are not None, are iterable,
        # and whose elements contain a "when" attribute.
        not_allowed = [
            "migmatrix",
            "asblack",
            "asdict",
            "fromdict",
            "_timed_events",
            "_state_asdict",
            "_reset_state",
        ]
        return [i for i in events.__dir__() if "__" not in i and i not in not_allowed]

    def _make_event_queues(self, events):
        """
        Take references from the input so that we can process them
        w/o affecting the input.
        """
        from collections import deque

        rv = {}
        for e in self._get_event_names(events):
            if events.__getattribute__(e) is not None:
                rv[e] = deque([i for i in events.__getattribute__(e)])
            else:
                rv[e] = deque([])

        return rv

    def _get_next_event_time(self, event_queues):
        t = None
        for _, value in event_queues.items():
            if len(value) > 0:
                if t is None:
                    t = value[0].when
                else:
                    t = min(t, value[0].when)
        return t

    def _current_events(self, t, event_queues, event):
        elist = []
        while len(event_queues[event]) > 0 and event_queues[event][0].when == t:
            e = event_queues[event][0]
            elist.append(e)
            event_queues[event].popleft()
        return elist

    def _label_deme(self, d):
        if self.deme_labels is None:
            return d
        if d not in self.deme_labels:
            return d
        return self.deme_labels[d]

    def _format_deme_sizes(self, sizes):
        rv = collections.OrderedDict()
        for i, j in enumerate(sizes):
            rv[self._label_deme(i)] = j
        return [(i, j) for i, j in rv.items()]

    def _apply_MassMigration(self, t, event_queues):
        for e in self._current_events(t, event_queues, "mass_migrations"):
            if self.current_deme_sizes[e.source] == 0:
                temp = (e, self._label_deme(e.source))
                raise ValueError(
                    "mass migration at time {} involves "
                    "empty source deme {}".format(*temp)
                )
            # NOTE: what is the C++ back-end doing if n_from_source < 1?
            # Answer: currently, that passes silently--is that okay?
            n_from_source = np.rint(
                self.current_deme_sizes[e.source] * e.fraction
            ).astype(int)
            if n_from_source == 0:
                # Assume that this occurs from rounding issues,
                # so we treat it as a warning rather than an error
                temp = (t, self._label_deme(e.source), self._label_deme(e.destination))
                message = (
                    "mass migration at time {} from {} " + "to {} moves no individuals"
                )
                warnings.warn(message.format(*temp))
            if (
                n_from_source > self.current_deme_sizes[e.source]
                and e.move_individuals is True
            ):
                temp = (t, n_from_source, self._label_deme(e.source))
                raise ValueError(
                    "mass migration at time {} "
                    "moves {} individuals from "
                    "source deme {}".format(*temp)
                )
            if e.move_individuals is True:
                self.current_deme_sizes[e.source] -= n_from_source
                self.current_deme_sizes[e.destination] += n_from_source
                temp = (
                    n_from_source,
                    self._label_deme(e.source),
                    self._label_deme(e.destination),
                )
                self._report.append(
                    "\tMass movement of {} " "from {} to {}\n".format(*temp)
                )
                if n_from_source > 0:
                    self.has_metadata[e.destination] = 1
                if e.fraction == 1.0:
                    self.went_extinct[e.source] = 1
                    self.has_metadata[e.source] = 0
            else:  # A copy
                if n_from_source > 0:
                    self.has_metadata[e.destination] = 1
                self.current_deme_sizes[e.destination] += n_from_source
                temp = (
                    n_from_source,
                    self._label_deme(e.source),
                    self._label_deme(e.destination),
                )
                self._report.append(
                    "\tMass copy of {} " "from {} to {}\n".format(*temp)
                )
            if e.resets_growth_rate is True:
                self.growth_rates[e.source] = fwdpy11.NOGROWTH
                self.growth_rates[e.destination] = fwdpy11.NOGROWTH
                s = "\t\tGrowth rates reset to {} in {} and {}\n"
                temp = (
                    fwdpy11.NOGROWTH,
                    self._label_deme(e.source),
                    self._label_deme(e.destination),
                )
                self._report.append(s.format(*temp))
            # Even if growth rates are not reset,
            # the onset times and initial sizes are affected
            self.growth_onset_times[e.source] = t
            self.growth_onset_times[e.destination] = t
            self.growth_initial_sizes[e.source] = self.current_deme_sizes[e.source]
            self.growth_initial_sizes[e.destination] = self.current_deme_sizes[
                e.destination
            ]
            self._report.append("\t\tGrowth onset times changed:\n")
            temp = (t, self._label_deme(e.source))
            self._report.append("\t\t\t{} in deme {}\n".format(*temp))
            temp = (t, self._label_deme(e.destination))
            self._report.append("\t\t\t{} in deme {}\n".format(*temp))
            self._report.append("\t\tGrowth initial sizes changed:\n")
            temp = (self.growth_initial_sizes[e.source], self._label_deme(e.source))
            self._report.append("\t\t\t{} in deme {}\n".format(*temp))
            temp = (
                self.growth_initial_sizes[e.destination],
                self._label_deme(e.destination),
            )
            self._report.append("\t\t\t{} in deme {}\n".format(*temp))

    def _apply_SetDemeSize(self, t, event_queues):
        for e in self._current_events(t, event_queues, "set_deme_sizes"):
            self.current_deme_sizes[e.deme] = e.new_size
            if e.new_size == 0:
                self.went_extinct[e.deme] = 1
            temp = (e.new_size, self._label_deme(e.deme))
            self._report.append("\tDeme size set to {} " "in deme {}\n".format(*temp))
            if e.resets_growth_rate is True:
                self.growth_rates[e.deme] = fwdpy11.NOGROWTH
                temp = (fwdpy11.NOGROWTH,)
                self._report.append("\t\tGrowth rate set to " "{}\n".format(*temp))

            # Deme size has change.  So, no matter what,
            # there is a new onsite time for growth!
            self.growth_onset_times[e.deme] = t
            self.growth_initial_sizes[e.deme] = self.current_deme_sizes[e.deme]
            temp = self.growth_initial_sizes[e.deme]
            self._report.append("\t\tGrowth initial size " "set to {}\n".format(temp))

    def _apply_SetSelfingRate(self, t, event_queues):
        for e in self._current_events(t, event_queues, "set_selfing_rates"):
            if self.current_deme_sizes[e.deme] == 0:
                temp = (self._label_deme(e.deme), e.when)
                raise ValueError(
                    "Setting selfing rate in "
                    "extinct deme {} at "
                    "time {}".format(*temp)
                )
            self.selfing_rates[e.deme] = e.S
            temp = (e.S, self._label_deme(e.deme))
            self._report.append(
                "\tSelfing probability " "set to {} in deme {}\n".format(*temp)
            )

    def _apply_SetExponentialGrowth(self, t, event_queues):
        for e in self._current_events(t, event_queues, "set_growth_rates"):
            if self.current_deme_sizes[e.deme] == 0:
                temp = (self._label_deme(e.deme), e.when)
                raise ValueError(
                    "attempt to change growth "
                    "rate in extinct deme {} "
                    "at time {}".format(*temp)
                )
            self.growth_rates[e.deme] = e.G
            self.growth_onset_times[e.deme] = t
            self.growth_initial_sizes[e.deme] = self.current_deme_sizes[e.deme]
            temp = (self._label_deme(e.deme), e.G)
            self._report.append("\tGrowth rate " "in deme {} set to {}\n".format(*temp))
            self._report.append("\t\tOnset time: {}\n".format(t))
            N = self.growth_initial_sizes[e.deme]
            self._report.append("\t\tOnset size: {}\n".format(N))

    def _apply_SetMigrationRates(self, t, event_queues):
        for e in self._current_events(t, event_queues, "set_migration_rates"):
            if e.deme >= 0:
                self._events.migmatrix.M[e.deme, :] = e.migrates
                temp = (self._label_deme(e.deme), e.migrates)
                self._report.append(
                    "\tMigration rates into " "deme {} set to {}\n".format(*temp)
                )
            else:
                self._events.migmatrix.M[:] = e.migrates.reshape(
                    self._events.migmatrix.M.shape
                )
                self._report.append(
                    "\tMigration matrix "
                    "reset to:\n\t\t{}".format(self._events.migmatrix.M)
                )

    def _apply_growth_rates(self, t, event_queues):
        next_deme_sizes = np.copy(self.current_deme_sizes)
        for i, j in enumerate(self.growth_rates):
            if j != fwdpy11.NOGROWTH:
                if self.current_deme_sizes[i] == 0:
                    temp = (t, self._label_deme(i))
                    raise ValueError(
                        "growth is happening "
                        "at time {} in extinct "
                        "deme {}".format(*temp)
                    )
                onset = self.growth_onset_times[i]
                G = self.growth_rates[i]
                Nonset = self.growth_initial_sizes[i]
                nextN = np.rint(Nonset * np.power(G, t - onset + 1)).astype(int)
                if nextN <= 0:
                    nextN = 0
                    self.growth_rates[i] = fwdpy11.NOGROWTH
                    self.growth_initial_sizes[i] = 0
                    self.growth_onset_times[i] = t
                    self.went_extinct[i] = 1
                next_deme_sizes[i] = nextN

        return next_deme_sizes

    def _validate_migration_rates(self, t, next_deme_sizes):
        if self._events.migmatrix is None:
            return
        for i, j in enumerate(next_deme_sizes):
            if j == 0:
                if self._events.migmatrix.M[i,].sum() != 0:
                    temp = (self._label_deme(i), t, self._events.migmatrix.M[i, j])
                    raise ValueError(
                        "there is migration "
                        "into empty deme {} "
                        "at time {}, "
                        "migrates={}".format(*temp)
                    )

    def _check_for_valid_parents(self, t):
        for i, j in enumerate(self.has_metadata):
            if j == 0 and self.current_deme_sizes[i] > 0:
                N = self.current_deme_sizes[i]
                if self._events.migmatrix is None:
                    s = "deme {} has size {} at time {} " + "but has no valid parents"
                    temp = (self._label_deme(i), N, t)
                    raise ValueError(s.format(*temp))
                elif self._events.migmatrix.M[i, i] != 0.0:
                    s = (
                        "deme {} at time {} "
                        + "has no valid parents in that deme "
                        + "but M[i, i] != 0"
                    )
                    temp = (self._label_deme(i), N, t)
                    raise ValueError(s.format(*temp))
                elif self._events.migmatrix.M[i,].sum() > 0:
                    self.has_metadata[i] = 1
                    self.went_extinct[i] = 0
            elif self.current_deme_sizes[i] > 0 and self._events.migmatrix is not None:
                # NOTE: Fix for issue 538, in 0.8.2
                if self._events.migmatrix.M[i,].sum() == 0:
                    s = (
                        f"deme {i} at time {t} has size {self.current_deme_sizes[i]} "
                        f"but an empty row in the migration matrix"
                    )
                    raise ValueError(s)

    def _generate_report(self, event_queues, simlen):
        """
        Apply events in the same way as the C++
        back-end.
        """
        temp = self._format_deme_sizes(self.current_deme_sizes)
        self._report = ["Deme sizes at time {}: {}\n".format(0, temp)]
        t = self._get_next_event_time(event_queues)
        if t > 0:
            # If the first event is after time zero, let's make
            # sure that the input migration matrix is valid.
            # This fixes github issue 544
            self._validate_migration_rates(0, self.current_deme_sizes)
        global_extinction = False
        last_t = None
        while t is not None and global_extinction is False:
            self.has_metadata[(self.went_extinct == 1)] = 0
            self.went_extinct[:] = 0
            self.current_deme_sizes[:] = self._apply_growth_rates(t - 1, event_queues)
            self._report.append("Events at time {}:\n".format(t))
            self._apply_MassMigration(t, event_queues)
            self._apply_SetDemeSize(t, event_queues)
            self._apply_SetExponentialGrowth(t, event_queues)
            self._apply_SetSelfingRate(t, event_queues)
            self._apply_SetMigrationRates(t, event_queues)
            next_deme_sizes = self._apply_growth_rates(t, event_queues)
            self._check_for_valid_parents(t)
            temp = self._format_deme_sizes(next_deme_sizes)
            sizes = "\tDeme sizes after growth: {}\n".format(temp)
            self._report.append(sizes)
            extinct = "\tThe following demes went extinct: {}\n"
            e = np.where(self.went_extinct == 1)[0]
            if len(e) > 0:
                self._report.append(extinct.format(e))
            extinct = "\tThe following demes are extinct: {}\n"
            e = np.where(self.has_metadata == 0)[0]
            if len(e) > 0:
                self._report.append(extinct.format(e))
            if next_deme_sizes.sum() == 0:
                global_extinction = True
                temp = "Global extinction occurs at time {}".format(t)
                warnings.warn(temp)
                self._report.append(temp + "\n")
            self._validate_migration_rates(t, next_deme_sizes)
            self.current_deme_sizes = next_deme_sizes
            last_t = t
            if simlen is not None and last_t > simlen:
                warnings.warn(f"current time {last_t} is > simulation length {simlen}")
            t = self._get_next_event_time(event_queues)

        if simlen is not None and simlen > last_t:
            deme_sizes = np.copy(self.current_deme_sizes)
            for i, j in enumerate(deme_sizes):
                if j > 0:
                    G = self.growth_rates[i]
                    if G != fwdpy11.NOGROWTH:
                        N0 = self.growth_initial_sizes[i]
                        t0 = self.growth_onset_times[i]
                        N = np.rint(N0 * np.power(G, simlen - t0)).astype(int)
                        deme_sizes[i] = N
            temp = "Final deme sizes at time {}: {}"
            deme_sizes = self._format_deme_sizes(deme_sizes)
            self._report.append(temp.format(simlen, deme_sizes))

    def _process_demographic_model(self, events, simlen):
        event_queues = self._make_event_queues(events)
        self._generate_report(event_queues, simlen)

    @property
    def report(self):
        """
        Obtain the details of the demographic
        model as a nicely-formatting string.
        """
        return str("").join(self._report)
