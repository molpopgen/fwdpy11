import math
from typing import Dict, Optional, Union

import attr
import demes
import numpy as np

from fwdpy11._types.demographic_model_citation import DemographicModelCitation
from fwdpy11._types.demographic_model_details import DemographicModelDetails
from fwdpy11._types.forward_demes_graph import ForwardDemesGraph


# TODO: need type hints for dg
def demography_from_demes(
    dg: Union[str, demes.Graph], burnin: int,
    round_non_integer_sizes=Optional[bool],
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

    fg = ForwardDemesGraph.from_demes(
        g, burnin, round_non_integer_sizes=round_non_integer_sizes)

    demography = _build_from_foward_demes_graph(fg, burnin, source)
    return demography


def _build_from_foward_demes_graph(
    fg: ForwardDemesGraph, burnin: int, source: Optional[Dict] = None
) -> DemographicModelDetails:
    """
    The workhorse.
    """
    idmap = _build_deme_id_to_int_map(fg.graph)
    initial_sizes = _get_initial_deme_sizes(fg.graph, idmap)
    Nref = _get_ancestral_population_size(fg.graph)

    burnin_generation = int(np.rint(burnin * Nref))
    model_times = _ModelTimes.from_demes_graph(fg.graph, burnin_generation)

    # TODO: size_history now contains model_times, so passing
    # the latter into functions is redundant.
    # We should clean this up later.
    # size_history = _DemeSizeHistory.from_demes_graph(
    #     dg, burnin, idmap, model_times)
    # assert size_history.model_times is not None

    # _set_initial_migration_matrix(dg, idmap, events, size_history)
    # _process_all_epochs(dg, idmap, model_times, events, size_history)
    # _process_migrations(dg, idmap, model_times, events, size_history)
    # _process_pulses(dg, idmap, model_times, events, size_history)
    # _process_admixtures(dg, dg_events, idmap,
    #                     model_times, events, size_history)
    # _process_mergers(dg, dg_events, idmap, model_times, events, size_history)
    # _process_splits(dg, dg_events, idmap, model_times, events, size_history)
    # _process_branches(dg, dg_events, idmap, model_times, events, size_history)

    if fg.graph.doi != "None":
        doi = fg.graph.doi
    else:
        doi = None

    return DemographicModelDetails(
        model=fg,
        name=fg.graph.description,
        source=source,
        parameters=None,
        citation=DemographicModelCitation(
            DOI=doi, full_citation=None, metadata=None),
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
    burnin_generation: int = attr.ib(
        validator=attr.validators.instance_of(int))

    @staticmethod
    def from_demes_graph(dg: demes.Graph, burnin_generation: int) -> "_ModelTimes":
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
            ends_inf = [
                d.epochs[0].end_time for d in dg.demes if d.start_time == math.inf
            ]
            starts = [d.start_time for d in dg.demes if d.start_time != math.inf]
            mig_starts = [
                m.start_time for m in dg.migrations if m.start_time != math.inf
            ]
            mig_ends = [
                m.end_time for m in dg.migrations if m.start_time == math.inf]
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
            burnin_generation=burnin_generation,
        )

    def convert_time(self, demes_event_time: float) -> int:
        """
        Backwards time -> forwards time
        """
        if demes_event_time != math.inf:
            return self.burnin_generation + int(
                self.model_start_time - demes_event_time - 1
            )

        return 0


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
    rate_change: float = attr.ib(
        validator=[attr.validators.instance_of(float)])
    from_deme_graph: bool = attr.ib(
        validator=attr.validators.instance_of(bool))


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
            rv[idmap[deme.name]] = int(np.rint(deme.epochs[0].start_size))

    if len(rv) == 0:
        raise RuntimeError("could not determine initial deme sizes")

    return rv


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

    temp = sorted(temp, key=lambda x: -x[0])

    return {j[1]: i for i, j in enumerate(temp)}


def _get_most_ancient_deme_start_time(dg: demes.Graph) -> demes.demes.Time:
    return max([d.start_time for d in dg.demes])


def _get_most_recent_deme_end_time(dg: demes.Graph) -> demes.demes.Time:
    return min([d.end_time for d in dg.demes])


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
            int(np.rint(e.start_size))
            for d in dg.demes
            for e in d.epochs
            if e.start_time == oldest_deme_time
        ]
    )
    if rv == 0:
        raise ValueError(
            "could not determinine ancestral metapopulation size")
    return rv
