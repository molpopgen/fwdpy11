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
from dataclasses import dataclass
import copy
import warnings
from typing import Dict, List, Optional, Union

import fwdpy11
import numpy as np

from ..demographic_models import DemographicModelDetails
from ..discrete_demography import (
    DiscreteDemography,
    _DemeSizeHistory,
    SetSelfingRate,
    SetDemeSize,
    SetExponentialGrowth,
    MassMigration,
    SetMigrationRates,
    MigrationMatrix,
)
from .diploid_population import DiploidPopulation


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
        sizes = o.deme_sizes(as_dict=True)
        return sizes
    except AttributeError:
        if isinstance(o, dict):
            return o
        else:
            rv = {}
            for i, s in enumerate(o):
                if s > 0:
                    rv[i] = s
            return rv


def _convert_simlen(simlen):
    if simlen is None:
        return simlen
    return int(simlen)


@dataclass
class _DemographyEvents:
    migmatrix: Optional[MigrationMatrix]
    mass_migrations: Optional[List[MassMigration]]
    set_deme_sizes: Optional[List[SetDemeSize]]
    set_growth_rates: Optional[List[SetExponentialGrowth]]
    set_selfing_rates: Optional[List[SetSelfingRate]]
    set_migration_rates: Optional[List[SetMigrationRates]]


@dataclass
class MigrationMatrixEpoch:
    start: int
    migration_matrix: np.ndarray


@dataclass
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

    initial_deme_sizes: Union[List[int], Dict[int, int], DiploidPopulation]
    events: Union[DiscreteDemography, DemographicModelDetails]
    simlen: Optional[int] = None
    deme_labels: Optional[Dict] = None

    def __post_init__(self):
        sizes = _create_initial_deme_sizes(self.initial_deme_sizes)
        if isinstance(self.events, DiscreteDemography):
            self._size_history = _DemeSizeHistory.from_lowlevel(
                sizes,
                mass_migrations=self.events.mass_migrations,
                set_deme_sizes=self.events.set_deme_sizes,
                set_growth_rates=self.events.set_growth_rates,
                total_simulation_length=self.simlen,
            )
            temp: fwdpy11.DiscreteDemography = self.events
        else:
            self._size_history = _DemeSizeHistory.from_lowlevel(
                sizes,
                mass_migrations=self.events.model.mass_migrations,
                set_deme_sizes=self.events.model.set_deme_sizes,
                set_growth_rates=self.events.model.set_growth_rates,
                total_simulation_length=self.simlen,
            )
            temp: fwdpy11.DiscreteDemography = self.events.model

        if temp.migmatrix is not None:
            event_list = _DemographyEvents(
                np.copy(temp.migmatrix.migmatrix),
                temp.mass_migrations,
                temp.set_deme_sizes,
                temp.set_growth_rates,
                temp.set_selfing_rates,
                temp.set_migration_rates,
            )
        else:
            event_list = _DemographyEvents(
                temp.migmatrix,
                temp.mass_migrations,
                temp.set_deme_sizes,
                temp.set_growth_rates,
                temp.set_selfing_rates,
                temp.set_migration_rates,
            )

        self.migration_history = self._build_migration_history(event_list)
        self._validate_epochs(event_list)

    def _build_migration_history(
        self, events: _DemographyEvents
    ):
        if events.migmatrix is None:
            return None

        current_migmatrix = np.copy(events.migmatrix)
        if current_migmatrix.ndim == 1 and len(current_migmatrix) == 1:
            warnings.warn(
                "You are using a 1x1 migration matrix."
                " You should prefer a migration matrix set to None."
            )

        # For the input migration matrix, make sure that all
        # source/dest pairings exist in generation 1
        for dest in range(current_migmatrix.shape[0]):
            for source, rate in enumerate(current_migmatrix[dest, :]):
                if rate > 0.0:
                    if self._size_history.deme_exists_at(source, 1) is False:
                        raise ValueError(
                            f"deme {source} does not exist at generation 1"
                        )
                    if self._size_history.deme_exists_at(dest, 1) is False:
                        raise ValueError(
                            f"deme {dest} does not exist at generation 1")

        rv = []
        if events.set_migration_rates is None:
            rv.append(MigrationMatrixEpoch(
                0, current_migmatrix))
        else:
            times = []
            for m in events.set_migration_rates:
                if m.when not in times:
                    times.append(m.when)
            times = sorted(times)
            for time in times:
                destinations = set()  # type: ignore
                for m in [m for m in events.set_migration_rates if m.when == time]:
                    if m.deme is not None and m.deme < 0:
                        if len(destinations) > 0:
                            raise ValueError(
                                f"entire migration matrix changed at when={m.when} in combination with"
                                " migration rate changes set for other demes individually"
                            )
                        if current_migmatrix.shape != m.migrates.shape:
                            raise ValueError(
                                f"migration matrix change at when={m.when} has an invalid shape"
                            )
                        current_migmatrix = np.copy(m.migrates).reshape(
                            current_migmatrix.shape
                        )
                    else:
                        if -1 in destinations:
                            raise ValueError(
                                f"entire migration matrix changed at when={m.when} in combination with"
                                " migration rate changes set for other demes individually"
                            )
                        if m.deme in destinations:
                            raise ValueError(
                                f"mutltiple migration rate changes for deme {m.deme} at when={m.when}"
                            )
                        # TODO: This will error out via numpy if the dimensions are bad.
                        # It would be better to have a clearer message.
                        current_migmatrix[m.deme, :] = np.copy(m.migrates)
                    destinations.add(m.deme)

                rv.append(MigrationMatrixEpoch(
                    time + 1, np.copy(current_migmatrix)))
        return rv

    def _validate_epochs(self, events: _DemographyEvents):
        for epoch in self._size_history.epochs:
            if epoch.begin > 0:
                when_ancestor_must_exist = epoch.begin - 1
            else:
                when_ancestor_must_exist = 0
            if epoch.data.ancestors is not None:
                for a in epoch.data.ancestors:
                    if (
                        self._size_history.deme_exists_at(
                            a, when_ancestor_must_exist)
                        is False
                    ):
                        raise ValueError(
                            f"deme {a} does not exist at {when_ancestor_must_exist}"
                        )
                    if (
                        when_ancestor_must_exist > 0
                        and self.migration_history is not None
                    ):
                        i = None
                        for j in reversed([m for m in self.migration_history]):
                            # TODO: is this correct?
                            # Go back and check how we are building
                            # the migration history
                            if when_ancestor_must_exist >= j.start:
                                i = j
                                break
                        if i is not None:
                            ttl_rate_in = i.migration_matrix[epoch.data.deme, :].sum(
                            )
                            if not ttl_rate_in > 0:
                                raise ValueError(
                                    f"the migration rantes into deme {epoch.data.deme} are 0.0 at time {i.start}"
                                )
            else:
                if events.migmatrix is None:
                    raise ValueError(
                        "a MigrationMatrix is required for this model.")
                else:
                    # TODO: does migration fix this?
                    if self.migration_history is not None:
                        i = None
                        for j in reversed([m for m in self.migration_history]):
                            if when_ancestor_must_exist >= j.start:
                                i = j
                                break
                        if i is not None:
                            row = i.migration_matrix[epoch.data.deme, :]
                            for source, rate in enumerate(row):
                                source_exists = self._size_history.deme_exists_at(
                                    source, i.start - 1
                                )
                                if rate > 0.0 and source_exists is False:
                                    raise ValueError(
                                        f"migration rate is > 0.0 from deme {source}, which doest not exist at time {i.start-1}"
                                    )

    @property
    def report(self):
        """
        Obtain the details of the demographic
        model as a nicely-formatting string.
        """
        warnings.warn("report generation has not been implemented")
        return "DemographyDebugger::report not yet implemented"
