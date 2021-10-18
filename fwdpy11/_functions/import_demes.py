import copy
import itertools
import math
import sys
from typing import Dict, List, Optional, Union

import attr
import demes
import numpy as np

from .. import class_decorators
from .._demography import exponential_growth_rate
from ..demographic_models import (DemographicModelCitation,
                                  DemographicModelDetails)
from ..discrete_demography import (DiscreteDemography, MassMigration,
                                   MigrationMatrix, SetDemeSize,
                                   SetExponentialGrowth, SetMigrationRates,
                                   SetSelfingRate, copy_individuals,
                                   move_individuals)


# TODO: need type hints for dg
def demography_from_demes(
    dg: Union[str, demes.Graph], burnin: int
) -> DemographicModelDetails:
    """
    The deme graph, dg, can be either a string or a resolved deme-graph.

    :param dg: A demes graph, either the object itself or a string for the YAML
        specified demography.
    :param int burnin: The factor for the burnin time, so that burn in occurs for
        burnin * N generations.

    .. versionadded:: 0.14.0
    """
    if isinstance(dg, str):
        g = demes.load(dg)
        source = {"demes_yaml_file": dg}
    else:
        g = dg
        source = None

    if not isinstance(burnin, int):
        raise ValueError("The burn in factor must be an integer")
    if burnin < 0:
        raise ValueError("Burn in factor must be non-negative")

    demography = _build_from_deme_graph(g, burnin, source)
    return demography


def _build_from_deme_graph(
    dg: demes.Graph, burnin: int, source: Optional[Dict] = None
) -> DemographicModelDetails:
    """
    The workhorse.
    """
    # dg must be in generations - replaces time_converter
    dg = dg.in_generations()

    # get classified discrete demographic events
    dg_events = dg.discrete_demographic_events()

    idmap = _build_deme_id_to_int_map(dg)
    initial_sizes = _get_initial_deme_sizes(dg, idmap)
    Nref = _get_ancestral_population_size(dg)

    burnin_generation = int(np.rint(burnin * Nref))
    model_times = _get_model_times(dg)

    events = _Fwdpy11Events(idmap=idmap)

    _set_initial_migration_matrix(dg, idmap, events)
    _process_all_epochs(dg, idmap, model_times, burnin_generation, events)
    _process_migrations(dg, idmap, model_times, burnin_generation, events)
    _process_pulses(dg, idmap, model_times, burnin_generation, events)
    _process_admixtures(dg, dg_events, idmap, model_times, burnin_generation, events)
    _process_mergers(dg, dg_events, idmap, model_times, burnin_generation, events)
    _process_splits(dg, dg_events, idmap, model_times, burnin_generation, events)
    _process_branches(dg, dg_events, idmap, model_times, burnin_generation, events)

    if dg.doi != "None":
        doi = dg.doi
    else:
        doi = None

    return DemographicModelDetails(
        model=events.build_model(),
        name=dg.description,
        source=source,
        parameters=None,
        citation=DemographicModelCitation(DOI=doi, full_citation=None, metadata=None),
        metadata={
            "deme_labels": {j: i for i, j in idmap.items()},
            "initial_sizes": initial_sizes,
            "burnin_time": burnin_generation,
            "total_simulation_length": burnin_generation
            + model_times.model_duration
            - 1,
        },
    )


@attr.s(frozen=True, auto_attribs=True)
class _ModelTimes(object):
    """
    These are in units of the deme graph
    and increase backwards into the past.
    """

    model_start_time: demes.demes.Time
    model_end_time: demes.demes.Time
    model_duration: int = attr.ib(validator=attr.validators.instance_of(int))


@attr.s(auto_attribs=True)
class _MigrationRateChange(object):
    """
    Use to make registry of migration rate changes.

    """

    when: int = attr.ib(
        validator=[demes.demes.non_negative, attr.validators.instance_of(int)]
    )
    source: int = attr.ib(
        validator=[demes.demes.non_negative, attr.validators.instance_of(int)]
    )
    destination: int = attr.ib(
        validator=[demes.demes.non_negative, attr.validators.instance_of(int)]
    )
    rate_change: float = attr.ib(validator=[attr.validators.instance_of(float)])
    from_deme_graph: bool = attr.ib(validator=attr.validators.instance_of(bool))


@class_decorators.attr_class_to_from_dict_no_recurse
@attr.s(auto_attribs=True)
class _Fwdpy11Events(object):
    """
    One stop shop for adding things we support in future versions.

    This class creates some redundancy with fwdpy11.DiscreteDemography,
    but that class is "frozen", so we need a mutable version.
    """

    mass_migrations: List[MassMigration] = attr.Factory(list)
    set_deme_sizes: List[SetDemeSize] = attr.Factory(list)
    set_growth_rates: List[SetExponentialGrowth] = attr.Factory(list)
    set_selfing_rates: List[SetSelfingRate] = attr.Factory(list)
    idmap: Dict = None

    # The initial continuous migration matrix
    initial_migmatrix: Optional[MigrationMatrix] = None
    # The migration matrix that we update to get changes in migration rates
    migmatrix: Optional[MigrationMatrix] = None

    # The following do not correspond to fwdpy11 event types.
    migration_rate_changes: List[_MigrationRateChange] = attr.Factory(list)
    # deme_extinctions: List[_DemeExtinctionEvent] = attr.Factory(list)

    def _update_changes_at_m(self, changes_at_m, migration_rate_change):
        # tally changes
        if migration_rate_change.from_deme_graph:
            changes_at_m[0][migration_rate_change.destination][
                migration_rate_change.source
            ] += migration_rate_change.rate_change
            if migration_rate_change.destination != migration_rate_change.source:
                changes_at_m[0][migration_rate_change.destination][
                    migration_rate_change.destination
                ] -= migration_rate_change.rate_change
        else:
            changes_at_m[1][migration_rate_change.destination][
                migration_rate_change.source
            ] += migration_rate_change.rate_change
            changes_at_m[1][migration_rate_change.destination][
                migration_rate_change.destination
            ] -= migration_rate_change.rate_change

    def _update_continuous_mass_migrations(self, changes_at_m, M_cont, M_mass):
        M_cont += changes_at_m[0]
        M_mass += changes_at_m[1]
        return M_cont, M_mass

    def _migration_matrix_from_partition(self, M_cont, M_mass):
        new_migmatrix = np.diag(np.diag(M_mass)).dot(M_cont) + (
            M_mass - np.diag(np.diag(M_mass))
        )
        return new_migmatrix

    def _build_migration_rate_changes(self) -> List[SetMigrationRates]:
        # We track the coninuous migration rates, and then augment with a matrix that
        # specifies changes to migration due to "instantaneous" events. The
        # instantaneous migration matrix is typically just the identity matrix,
        # (i.e. uses continuous rates, but off diag elements scale continuous rates
        # and add to ancestry source from the off diagonal source column)
        # but for some generations has ancestry pointing to different demes due
        # to pulse, split, etc events
        if self.migmatrix is not None:
            M_cont = copy.deepcopy(self.migmatrix)
        else:
            M_cont = np.eye(len(self.idmap))
        M_mass = np.eye(len(self.idmap))

        set_migration_rates: List[SetMigrationRates] = []

        self.migration_rate_changes = sorted(
            self.migration_rate_changes,
            key=lambda x: (x.when, x.destination, x.source, x.from_deme_graph),
        )
        # self.deme_extinctions = sorted(
        #    self.deme_extinctions, key=lambda x: (x.when, x.deme)
        # )
        m = 0
        changes_at_m = [
            np.zeros((len(self.idmap), len(self.idmap))),  # from Graph
            np.zeros((len(self.idmap), len(self.idmap))),  # not from Graph
        ]
        while m < len(self.migration_rate_changes):
            # gather all migration rate changes and extinction events
            # that occur at a given time
            self._update_changes_at_m(changes_at_m, self.migration_rate_changes[m])
            mm = m + 1

            while (
                mm < len(self.migration_rate_changes)
                and self.migration_rate_changes[mm].when
                == self.migration_rate_changes[m].when
            ):
                self._update_changes_at_m(changes_at_m, self.migration_rate_changes[mm])
                mm += 1

            # update M_cont and M_mass
            M_cont, M_mass = self._update_continuous_mass_migrations(
                changes_at_m, M_cont, M_mass
            )
            # get the new migration matrix
            new_migmatrix = self._migration_matrix_from_partition(M_cont, M_mass)
            # for any rows that don't match, add a fwdpy11.SetMigrationRate
            for i in range(len(self.idmap)):
                if np.any(self.migmatrix[i] != new_migmatrix[i]):
                    set_migration_rates.append(
                        SetMigrationRates(
                            self.migration_rate_changes[m].when, i, new_migmatrix[i]
                        )
                    )

            self.migmatrix = new_migmatrix

            m = mm
            # reset changes
            changes_at_m = [
                np.zeros((len(self.idmap), len(self.idmap))),  # from Graph
                np.zeros((len(self.idmap), len(self.idmap))),  # not from Graph
            ]

        return set_migration_rates

    ### now unused ??
    def build_model(self) -> DiscreteDemography:
        set_migration_rates = self._build_migration_rate_changes()
        return DiscreteDemography(
            mass_migrations=self.mass_migrations,
            set_deme_sizes=self.set_deme_sizes,
            set_growth_rates=self.set_growth_rates,
            set_selfing_rates=self.set_selfing_rates,
            migmatrix=self.initial_migmatrix,
            set_migration_rates=set_migration_rates,
        )


def _build_deme_id_to_int_map(dg: demes.Graph) -> Dict:
    """
    Convert the string input ID to output integer values.

    For sanity, the output values will be in increasing
    order of "deme origin" times.

    We rely on the epoch times being strictly sorted
    past-to-present in the YAML.  If there are ties,
    populations are sorted lexically by ID within
    start_time.
    """
    temp = []
    for deme in dg.demes:
        assert len(deme.epochs) > 0
        temp.append((deme.epochs[0].start_time, deme.name))

    temp = sorted(temp, key=lambda x: (-x[0], x[1]))

    return {j[1]: i for i, j in enumerate(temp)}


def _get_initial_deme_sizes(dg: demes.Graph, idmap: Dict) -> Dict:
    """
    Build a map of a deme's integer label to its size
    at the start of the simulation for all demes whose
    start_time equals inf.
    """
    otime = _get_most_ancient_deme_start_time(dg)
    rv = dict()
    for deme in dg.demes:
        if deme.epochs[0].start_time == otime:
            rv[idmap[deme.name]] = deme.epochs[0].start_size

    if len(rv) == 0:
        raise RuntimeError("could not determine initial deme sizes")

    return rv


def _get_most_ancient_deme_start_time(dg: demes.Graph) -> demes.demes.Time:
    return max([d.start_time for d in dg.demes])


def _get_most_recent_deme_end_time(dg: demes.Graph) -> demes.demes.Time:
    return min([d.end_time for d in dg.demes])


def _get_model_times(dg: demes.Graph) -> _ModelTimes:
    """
    In units of dg.time_units, obtain the following:

    1. The time when the demographic model starts.
    2. The time when it ends.
    3. The total simulation length.

    """
    # FIXME: this function isn't working well.
    # For example, twodemes.yml and twodemes_one_goes_away.yml
    # both break it.
    oldest_deme_time = _get_most_ancient_deme_start_time(dg)
    most_recent_deme_end = _get_most_recent_deme_end_time(dg)

    model_start_time = oldest_deme_time
    if oldest_deme_time == math.inf:
        # We want to find the time of first event or
        # the first demographic change, which is when
        # burnin will end. To do this, get a list of
        # first size change for all demes with inf
        # start time, and the start time for all other
        # demes, and take max of those.
        ends_inf = [d.epochs[0].end_time for d in dg.demes if d.start_time == math.inf]
        starts = [d.start_time for d in dg.demes if d.start_time != math.inf]
        mig_starts = [m.start_time for m in dg.migrations if m.start_time != math.inf]
        mig_ends = [m.end_time for m in dg.migrations if m.start_time == math.inf]
        pulse_times = [p.time for p in dg.pulses]
        # The forward-time model with start with a generation 0,
        # which is the earliest end point of a deme with start time
        # of inf, minus 1.  That definition is forwards in time, so we
        # ADD one to the backwards-in-time demes info.
        model_start_time = (
            max(ends_inf + starts + mig_starts + mig_ends + pulse_times) + 1
        )

    if most_recent_deme_end != 0:
        model_duration = model_start_time - most_recent_deme_end
    else:
        model_duration = model_start_time

    return _ModelTimes(
        model_start_time=model_start_time,
        model_end_time=most_recent_deme_end,
        model_duration=int(np.rint(model_duration)),
    )


def _get_ancestral_population_size(dg: demes.Graph) -> int:
    """
    Need this for the burnin time.

    If there are > 1 demes with the same most ancient start_time,
    then the ancestral size is considered to be the size
    of all those demes (size of ancestral metapopulation).
    """
    oldest_deme_time = _get_most_ancient_deme_start_time(dg)

    rv = sum(
        [
            e.start_size
            for d in dg.demes
            for e in d.epochs
            if e.start_time == oldest_deme_time
        ]
    )
    if rv == 0:
        raise RuntimeError("could not determinine ancestral metapopulation size")
    return rv


def _set_initial_migration_matrix(
    dg: demes.Graph, idmap: Dict, events: _Fwdpy11Events
) -> None:
    """
    Set any migration rates that have start time of inf. More
    recent migration rate changes or start time set the model
    start time at or before that time
    """
    if len(idmap) > 1:
        migmatrix = np.zeros((len(idmap), len(idmap)))
        for deme_id, ii in idmap.items():
            if dg[deme_id].start_time == math.inf:
                migmatrix[ii, ii] = 1.0
        if len(dg.migrations) > 0:
            for m in dg.migrations:
                if m.start_time == math.inf:
                    migmatrix[idmap[m.dest]][idmap[m.source]] += m.rate
                    migmatrix[idmap[m.dest]][idmap[m.dest]] -= m.rate

        events.migmatrix = migmatrix
        events.initial_migmatrix = migmatrix


##
## Process deme epochs, migrations, and discrete demographic events
##


def _process_epoch(
    deme_id: str,
    e: demes.Epoch,
    idmap: Dict,
    model_times: _ModelTimes,
    burnin_generation: int,
    events: _Fwdpy11Events,
) -> None:
    """
    Can change sizes, cloning rates, and selfing rates.

    Since fwdpy11 currently doesn't understand cloning, we need
    to raise an error if the rate is not None or nonzero.
    """
    if e.start_time != math.inf:
        when = burnin_generation + int(model_times.model_start_time - e.start_time - 1)
    else:
        when = 0

    if e.selfing_rate is not None:
        events.set_selfing_rates.append(
            SetSelfingRate(when=when, deme=idmap[deme_id], S=e.selfing_rate)
        )

    if e.cloning_rate is not None:
        if e.cloning_rate > 0.0:
            raise ValueError("fwdpy11 does not currently support cloning rates > 0.")

    if e.start_time != math.inf:
        # Handle size change functions
        events.set_deme_sizes.append(
            SetDemeSize(
                when=when,
                deme=idmap[deme_id],
                new_size=e.start_size,
            )
        )
        if e.end_size != e.start_size:
            if e.size_function != "exponential":
                raise ValueError(
                    f"Size change function must be exponential.  We got {e.size_function}"
                )
            G = exponential_growth_rate(
                e.start_size, e.end_size, int(np.rint(e.time_span))
            )
            events.set_growth_rates.append(
                SetExponentialGrowth(when=when, deme=idmap[deme_id], G=G)
            )


def _process_all_epochs(
    dg: demes.Graph,
    idmap: Dict,
    model_times: _ModelTimes,
    burnin_generation: int,
    events: _Fwdpy11Events,
):
    """
    Processes all epochs of all demes to set sizes and selfing rates.
    """
    for deme in dg.demes:
        for e in deme.epochs:
            _process_epoch(
                deme.name,
                e,
                idmap,
                model_times,
                burnin_generation,
                events,
            )
        # if a deme starts more recently than math.inf, we have to
        # turn on migration in that deme with a diagonal element to 1
        if deme.start_time < math.inf:
            events.migration_rate_changes.append(
                _MigrationRateChange(
                    when=burnin_generation
                    + int(model_times.model_start_time - deme.start_time - 1),
                    source=idmap[deme.name],
                    destination=idmap[deme.name],
                    rate_change=1.0,
                    from_deme_graph=True,
                )
            )

        # if a deme ends before time zero, we set diag entry in migmatrix to 0
        if deme.end_time > 0:
            events.migration_rate_changes.append(
                _MigrationRateChange(
                    when=burnin_generation
                    + int(model_times.model_start_time - deme.end_time - 1),
                    source=idmap[deme.name],
                    destination=idmap[deme.name],
                    rate_change=-1.0,
                    from_deme_graph=True,
                )
            )

        # if deme ends before time zero, we set set its size to zero
        # we proces deme extinctions here instead of in the events
        if deme.end_time > 0:
            events.set_deme_sizes.append(
                SetDemeSize(
                    when=burnin_generation
                    + int(model_times.model_start_time - deme.end_time - 1),
                    deme=idmap[deme.name],
                    new_size=0,
                )
            )

        # collect all (deme, time) tuples that represent extinctions
        extinctions = [
            (i.deme, i.when) for i in events.set_deme_sizes if i.new_size == 0
        ]

        # purge invalid events:
        # * setting selfing rates in extinct demes
        # * setting growth rates in extinct demes
        events.set_selfing_rates = [
            i for i in events.set_selfing_rates if (i.deme, i.when) not in extinctions
        ]
        events.set_growth_rates = [
            i for i in events.set_growth_rates if (i.deme, i.when) not in extinctions
        ]


def _process_migrations(
    dg: demes.Graph,
    idmap: Dict,
    model_times: _ModelTimes,
    burnin_generation: int,
    events: _Fwdpy11Events,
) -> None:
    """
    Make a record of everything in dg.migrations

    When a migration rate has an end time > 0, it gets entered twice.
    """
    for m in dg.migrations:
        if m.start_time < math.inf:
            when = burnin_generation + int(
                model_times.model_start_time - m.start_time - 1
            )
            try:
                events.migration_rate_changes.append(
                    _MigrationRateChange(
                        when=when,
                        source=idmap[m.source],
                        destination=idmap[m.dest],
                        rate_change=m.rate,
                        from_deme_graph=True,
                    )
                )
            except AttributeError:
                for source, dest in itertools.permutations(m.demes, 2):
                    events.migration_rate_changes.append(
                        _MigrationRateChange(
                            when=when,
                            source=idmap[source],
                            destination=idmap[dest],
                            rate_change=m.rate,
                            from_deme_graph=True,
                        )
                    )
        if m.end_time > 0:
            when = burnin_generation + int(
                model_times.model_start_time - m.end_time - 1
            )
            try:
                events.migration_rate_changes.append(
                    _MigrationRateChange(
                        when=when,
                        source=idmap[m.source],
                        destination=idmap[m.dest],
                        rate_change=-m.rate,
                        from_deme_graph=True,
                    )
                )
            except AttributeError:
                for source, dest in itertools.permutations(m.demes, 2):
                    events.migration_rate_changes.append(
                        _MigrationRateChange(
                            when=when,
                            source=idmap[source],
                            destination=idmap[dest],
                            rate_change=-m.rate,
                            from_deme_graph=True,
                        )
                    )


def _process_pulses(
    dg: demes.Graph,
    idmap: Dict,
    model_times: _ModelTimes,
    burnin_generation: int,
    events: _Fwdpy11Events,
) -> None:
    for p in dg.pulses:
        when = burnin_generation + int(model_times.model_start_time - p.time - 1)
        events.migration_rate_changes.append(
            _MigrationRateChange(
                when=when,
                source=idmap[p.source],
                destination=idmap[p.dest],
                rate_change=p.proportion,
                from_deme_graph=False,
            )
        )
        events.migration_rate_changes.append(
            _MigrationRateChange(
                when=when + 1,
                source=idmap[p.source],
                destination=idmap[p.dest],
                rate_change=-p.proportion,
                from_deme_graph=False,
            )
        )


def _process_admixtures(
    dg: demes.Graph,
    dg_events: Dict,
    idmap: Dict,
    model_times: _ModelTimes,
    burnin_generation: int,
    events: _Fwdpy11Events,
) -> None:
    for a in dg_events["admixtures"]:
        when = burnin_generation + int(model_times.model_start_time - a.time - 1)
        for parent, proportion in zip(a.parents, a.proportions):
            events.migration_rate_changes.append(
                _MigrationRateChange(
                    when=when,
                    source=idmap[parent],
                    destination=idmap[a.child],
                    rate_change=proportion,
                    from_deme_graph=False,
                )
            )
            events.migration_rate_changes.append(
                _MigrationRateChange(
                    when=when + 1,
                    source=idmap[parent],
                    destination=idmap[a.child],
                    rate_change=-proportion,
                    from_deme_graph=False,
                )
            )


def _process_mergers(
    dg: demes.Graph,
    dg_events: Dict,
    idmap: Dict,
    model_times: _ModelTimes,
    burnin_generation: int,
    events: _Fwdpy11Events,
) -> None:
    for m in dg_events["mergers"]:
        when = burnin_generation + int(model_times.model_start_time - m.time - 1)
        for parent, proportion in zip(m.parents, m.proportions):
            events.migration_rate_changes.append(
                _MigrationRateChange(
                    when=when,
                    source=idmap[parent],
                    destination=idmap[m.child],
                    rate_change=proportion,
                    from_deme_graph=False,
                )
            )
            events.migration_rate_changes.append(
                _MigrationRateChange(
                    when=when + 1,
                    source=idmap[parent],
                    destination=idmap[m.child],
                    rate_change=-proportion,
                    from_deme_graph=False,
                )
            )


def _process_splits(
    dg: demes.Graph,
    dg_events: Dict,
    idmap: Dict,
    model_times: _ModelTimes,
    burnin_generation: int,
    events: _Fwdpy11Events,
) -> None:
    """
    A split is a "sudden" creation of > 1 offspring deme
    from a parent and the parent ceases to exist.

    Given that there is no proportions attribute, we infer/assume
    (danger!) that each offspring deme gets 100% of its ancestry
    from the parent.
    """
    for s in dg_events["splits"]:
        when = burnin_generation + int(model_times.model_start_time - s.time - 1)
        for c in s.children:
            # one generation of migration to move lineages from parent to children
            events.migration_rate_changes.append(
                _MigrationRateChange(
                    when=when,
                    source=idmap[s.parent],
                    destination=idmap[c],
                    rate_change=1.0,
                    from_deme_graph=False,
                )
            )
            # turn off that migration after one generation
            events.migration_rate_changes.append(
                _MigrationRateChange(
                    when=when + 1,
                    source=idmap[s.parent],
                    destination=idmap[c],
                    rate_change=-1.0,
                    from_deme_graph=False,
                )
            )


def _process_branches(
    dg: demes.Graph,
    dg_events: Dict,
    idmap: Dict,
    model_times: _ModelTimes,
    burnin_generation: int,
    events: _Fwdpy11Events,
) -> None:
    """
    A branch creates a child deme with 100% ancestry from the parent.
    The parent continues to exist.

    The 1-to-1 relationship between parent and child means 100% of the
    child's ancestry is from parent.
    """
    for b in dg_events["branches"]:
        when = burnin_generation + int(model_times.model_start_time - b.time - 1)
        # turn on migration for one generation at "when"
        events.migration_rate_changes.append(
            _MigrationRateChange(
                when=when,
                source=idmap[b.parent],
                destination=idmap[b.child],
                rate_change=1.0,
                from_deme_graph=False,
            )
        )
        # end that migration after one generation
        events.migration_rate_changes.append(
            _MigrationRateChange(
                when=when + 1,
                source=idmap[b.parent],
                destination=idmap[b.child],
                rate_change=-1.0,
                from_deme_graph=False,
            )
        )
