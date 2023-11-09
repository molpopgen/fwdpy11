#
# Copyright(C) 2017 Kevin Thornton < krthornt @uci.edu>
#
# This file is part of fwdpy11.
#
# fwdpy11 is free software : you can redistribute it and / or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# fwdpy11 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with fwdpy11.If not, see < http: //www.gnu.org/licenses/>.
#

from typing import Callable, Optional
import warnings

import fwdpy11

from ._fwdpy11 import GSLrng, SampleRecorder
from ._types import DiploidPopulation, ModelParams


def evolvets(
    rng: GSLrng,
    pop: DiploidPopulation,
    params: ModelParams,
    simplification_interval: int,
    recorder: Optional[Callable[[
        DiploidPopulation, SampleRecorder], None]] = None,
    *,
    post_simplification_recorder: Optional[Callable[[
        DiploidPopulation], None]] = None,
    suppress_table_indexing: Optional[bool] = None,
    record_gvalue_matrix: bool = False,
    stopping_criterion: Optional[Callable[[
        DiploidPopulation, bool], bool]] = None,
    track_mutation_counts: bool = False,
    remove_extinct_variants: bool = True,
    preserve_first_generation: bool = False,
):
    """
    Evolve a population with tree sequence recording

    :param rng: random number generator
    :type rng: :class:`fwdpy11.GSLrng`
    :param pop: A population
    :type pop: :class:`fwdpy11.DiploidPopulation`
    :param params: simulation parameters
    :type params: :class:`fwdpy11.ModelParams`
    :param simplification_interval: Number of generations
        between simplifications.
    :type simplification_interval: int
    :param recorder: (None) A temporal sampler/data recorder.
    :type recorder: typing.Callable
    :param post_simplification_recorder: (None) A temporal sampler
    :type post_simplification_recorder: typing.Callable
    :param suppress_table_indexing: (`None`) Prevents edge table indexing until
                                    end of simulation.  The default value
                                    (`None`) will be interpreted as `True`
    :type suppress_table_indexing: Optional[bool]
    :param record_gvalue_matrix: (False) Whether to record genetic values into
                                 :attr:`fwdpy11.DiploidPopulation.genetic_values`.
    :type record_gvalue_matrix: bool
    :param preserve_first_generation: (False) Whether to record generation 0 as
                                      ancient samples. Must be `True` for
                                      tree sequence "recapitation". See
                                      :ref:`recapitation`.
    :type preserve_first_generation: bool

    The recording of genetic values into :attr:`fwdpy11.DiploidPopulation.genetic_values`
    is suppressed by default.  First, it is redundant with
    :attr:`fwdpy11.DiploidMetadata.g` for the common case of mutational effects on a
    single trait. Second, we save some memory by not tracking these matrices.
    However, it is useful to track these data for some
    cases when simulating multivariate mutational effects (pleiotropy).

    .. note::
        If recorder is None,
        then :class:`fwdpy11.NoAncientSamples` will be used.

        If post_simplification_recorder is None, then
        :class:`fwdpy11.RecordNothing` will be used.

    .. versionchanged:: 0.5.2

        Added post_simplification_recorder.

    .. versionchanged:: 0.7.1

        Added preserve_first_generation.

    .. versionchanged:: 0.8.0

        Update to refactored ModelParams.
        Added ``check_demographic_event_timings``.

    .. versionchanged: 0.18.0

        `suppress_table_indexing` default changed to `None`,
        which is interpreted as `True`.

    .. versionchanged: 0.20.0

        * Back-end refactored to use demes models directly.
        * Will define a default demographic model using
          :func:`fwdpy11.DemesForwardGraph.tubes` if the
          input model is `None`.
        * Remove option `check_demographic_event_timings`

    """
    if params.demography is not None:
        try:
            # DemographicModelDetails ?
            demographic_model = params.demography.model
        except AttributeError:
            # No, assume it is a DiscreteDemography
            demographic_model = params.demography
        except Exception as e:
            raise e
    else:
        # Build a default model of "tubes"
        from fwdpy11 import ForwardDemesGraph
        sizes = pop.deme_sizes()[1].tolist()
        msg = "Applying a default demographic model "
        msg += f"where deme sizes are {sizes} "
        msg += f"and the burn-in length is 10*{sum(sizes)}. "
        msg += "This will raise an error in future releases."
        warnings.warn(msg, DeprecationWarning, stacklevel=2)
        demographic_model = ForwardDemesGraph.tubes(sizes, 10)

    for r in params.sregions:
        if isinstance(r, fwdpy11.mvDES):
            if isinstance(r.des, list):
                nd = demographic_model.number_of_demes()
                if nd > params.gvalue.ndemes:
                    msg = f"maxmimum number of demes in model ({nd})"
                    msg += " is not compatible with genetic value"
                    msg += f" number of demes {params.gvalue.ndemes}"
                    raise ValueError(msg)
                for rr in r.des:
                    assert isinstance(rr, fwdpy11.Sregion)  # type: ignore
                    if rr.beg < 0:
                        raise ValueError(f"{r} has begin value < 0.0")
                    if rr.end > pop.tables.genome_length:
                        raise ValueError(
                            f"{r} has end value >= genome length"
                            " of {pop.tables.genome_length}"
                        )
            else:
                assert isinstance(r, fwdpy11.mvDES), f"{type(r)}"
                assert isinstance(r.des, fwdpy11.Sregion)  # type: ignore
                nd = demographic_model.number_of_demes()
                if isinstance(params.gvalue, list):
                    for g in params.gvalue:
                        if nd > g.ndemes:
                            msg = f"maxmimum number of demes in model ({nd})"
                            msg += " is not compatible with genetic value"
                            msg += f" number of demes {params.gvalue.ndemes}"
                            raise ValueError(msg)
                else:
                    if nd > params.gvalue.ndemes:
                        msg = f"maxmimum number of demes in model ({nd})"
                        msg += " is not compatible with genetic value"
                        msg += f" number of demes {params.gvalue.ndemes}"
                        raise ValueError(msg)
                if r.des.beg < 0:
                    raise ValueError(f"{r} has begin value < 0.0")
                if r.des.end > pop.tables.genome_length:
                    raise ValueError(
                        f"{r} has end value >= genome"
                        " length of {pop.tables.genome_length}"
                    )
        else:
            assert isinstance(r, fwdpy11.Sregion)  # type: ignore
            if r.beg < 0:
                raise ValueError(f"{r} has begin value < 0.0")
            if r.end > pop.tables.genome_length:
                raise ValueError(
                    f"{r} has end value >= genome"
                    " length of {pop.tables.genome_length}"
                )

    for r in params.nregions:
        assert isinstance(r, fwdpy11.Region)
        if r.beg < 0:
            raise ValueError(f"{r} has begin value < 0.0")
        if r.end > pop.tables.genome_length:
            raise ValueError(
                f"{r} has end value > genome length of {pop.tables.genome_length}"
            )

    for r in params.recregions:
        assert isinstance(r, fwdpy11._fwdpy11.PoissonCrossoverGenerator) or \
            isinstance(
                r, fwdpy11._fwdpy11.NonPoissonCrossoverGenerator) or \
            isinstance(r, fwdpy11.Region), f"{type(r)}"
        if hasattr(r, "beg") and hasattr(r, "end"):
            if r.beg < 0:
                raise ValueError(f"{r} has begin value < 0.0")
            if r.end > pop.tables.genome_length:
                raise ValueError(
                    f"{r} has end value >= genome length of {pop.tables.genome_length}"
                )
        elif hasattr(r, "position"):
            if r.position < 0:
                raise ValueError(f"{r} has position value < 0.0")
            if r.position >= pop.tables.genome_length:
                raise ValueError(
                    f"{r} has position value >= genome length of {pop.tables.genome_length}"
                )
        elif hasattr(r, "regions"):
            for i in r.regions:
                if i.beg < 0:
                    raise ValueError(f"{r} has begin value < 0.0")
                if i.end > pop.tables.genome_length:
                    raise ValueError(
                        f"{r} has end value >= genome length of {pop.tables.genome_length}"
                    )
        else:
            raise TypeError(f"unexpceted type: {type(r)}")

    if recorder is None:
        from ._fwdpy11 import NoAncientSamples

        recorder = NoAncientSamples()

    if post_simplification_recorder is None:
        from ._fwdpy11 import RecordNothing

        post_simplification_recorder = RecordNothing()
        reset_treeseqs_after_simplify = False
    else:
        reset_treeseqs_after_simplify = True

    if stopping_criterion is None:
        from ._fwdpy11 import _no_stopping

        stopping_criterion = _no_stopping

    from ._fwdpy11 import (
        MutationRegions,
        dispatch_create_GeneticMap,
        dispatch_create_GeneticMap_non_Region,
        evolve_with_tree_sequences,
    )

    fwdpy11._validate_regions(params.sregions, pop.tables.genome_length)
    fwdpy11._validate_regions(params.nregions, pop.tables.genome_length)

    pneutral = 0.0
    if params.rates.neutral_mutation_rate + params.rates.selected_mutation_rate > 0.0:
        pneutral = params.rates.neutral_mutation_rate / (
            params.rates.neutral_mutation_rate + params.rates.selected_mutation_rate
        )
    mm = MutationRegions.create(pneutral, params.nregions, params.sregions)

    if all([isinstance(i, fwdpy11.Region) for i in params.recregions]):
        rm = dispatch_create_GeneticMap(
            params.rates.recombination_rate, params.recregions)
    else:
        poisson = [i for i in params.recregions if isinstance(
            i, fwdpy11._fwdpy11.PoissonCrossoverGenerator)]
        non_poisson = [i for i in params.recregions if isinstance(
            i, fwdpy11._fwdpy11.NonPoissonCrossoverGenerator)]
        rm = dispatch_create_GeneticMap_non_Region(poisson, non_poisson)
        assert rm._num_poisson_callbacks() == len(poisson)
        assert rm._num_non_poisson_callbacks() == len(non_poisson)

    from ._fwdpy11 import SampleRecorder, _evolve_with_tree_sequences_options

    sr = SampleRecorder()
    from ._fwdpy11 import _dgvalue_pointer_vector

    options = _evolve_with_tree_sequences_options()

    options.preserve_selected_fixations = params.prune_selected is False
    if suppress_table_indexing is not None:
        options.suppress_edge_table_indexing = suppress_table_indexing
    else:
        options.suppress_edge_table_indexing = True
    options.record_gvalue_matrix = record_gvalue_matrix
    options.track_mutation_counts_during_sim = track_mutation_counts
    options.remove_extinct_mutations_at_finish = remove_extinct_variants
    options.reset_treeseqs_to_alive_nodes_after_simplification = \
        reset_treeseqs_after_simplify
    options.preserve_first_generation = preserve_first_generation
    options.allow_residual_selfing = params.allow_residual_selfing

    if options.allow_residual_selfing is False:
        from fwdpy11._types.forward_demes_graph import _round_via_decimal
        for d, deme in enumerate(demographic_model.graph.demes):
            for e, epoch in enumerate(deme.epochs):
                start_size = _round_via_decimal(epoch.start_size)
                end_size = _round_via_decimal(epoch.end_size)
                if start_size == 1 or end_size == 1:
                    msg = f"deme {d} epoch {e} has a size of 1"
                    msg += " but residual selfing is not allowed."
                    msg += " This model will raise an error if this deme "
                    msg += "contributes to ancestry after the size of 1"
                    msg += " is reached."
                    warnings.warn(msg, UserWarning, stacklevel=2)

    gvpointers = _dgvalue_pointer_vector(params.gvalue)
    evolve_with_tree_sequences(
        rng,
        pop,
        sr,
        simplification_interval,
        demographic_model,
        params.simlen,
        params.rates.neutral_mutation_rate,
        params.rates.selected_mutation_rate,
        mm,
        rm,
        gvpointers,
        recorder,
        stopping_criterion,
        post_simplification_recorder,
        options,
    )
