from typing import Dict, Optional, Union

import demes

from fwdpy11._types.demographic_model_citation import DemographicModelCitation
from fwdpy11._types.demographic_model_details import DemographicModelDetails
from fwdpy11._types.forward_demes_graph import ForwardDemesGraph


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

    if fg.graph.doi != "None":
        doi = fg.graph.doi
    else:
        doi = None

    _initial_sizes = {i: j for i, j in enumerate(
        fg._parental_deme_sizes_at_time_zero()) if j > 0}
    return DemographicModelDetails(
        model=fg,
        name=fg.graph.description,
        source=source,
        parameters=None,
        citation=DemographicModelCitation(
            DOI=doi, full_citation=None, metadata=None),
        metadata={
            "deme_labels": {j: i for i, j in idmap.items()},
            "initial_sizes": _initial_sizes,
            "burnin_time": fg.burnin_generation,
            "total_simulation_length": fg._model_end_time() - 1
        },
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

    temp = sorted(temp, key=lambda x: -x[0])

    return {j[1]: i for i, j in enumerate(temp)}
