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

import fwdpy11


def _validate_event_timings(demography: fwdpy11.DiscreteDemography, generation: int):
    too_early = []
    for i in demography._timed_events():
        if i is not None:
            for j in i:
                if j.when < generation:
                    too_early.append(j)
    if len(too_early) > 0:
        import warnings

        msg = "The following demographic events are registered to occur "
        msg += f"before the current generation, which is {generation}: {too_early}"
        warnings.warn(msg)


def evolvets(
    rng,
    pop,
    params,
    simplification_interval,
    recorder=None,
    *,
    post_simplification_recorder=None,
    suppress_table_indexing=False,
    record_gvalue_matrix=False,
    stopping_criterion=None,
    track_mutation_counts=False,
    remove_extinct_variants=True,
    preserve_first_generation=False,
    check_demographic_event_timings=True,
):
    """
    Evolve a population with tree sequence recording

    :param rng: random number generator
    :type rng: :class:`fwdpy11.GSLrng`
    :param pop: A population
    :type pop: :class:`fwdpy11.DiploidPopulation`
    :param params: simulation parameters
    :type params: :class:`fwdpy11.ModelParams`
    :param simplification_interval: Number of generations between simplifications.
    :type simplification_interval: int
    :param recorder: (None) A temporal sampler/data recorder.
    :type recorder: typing.Callable
    :param post_simplification_recorder: (None) A temporal sampler
    :type post_simplification_recorder: typing.Callable
    :param suppress_table_indexing: (False) Prevents edge table indexing until
                                    end of simulation
    :type suppress_table_indexing: bool
    :param record_gvalue_matrix: (False) Whether to record genetic values into
                                 :attr:`fwdpy11.PopulationBase.genetic_values`.
    :type record_gvalue_matrix: bool
    :param preserve_first_generation: (False) Whether to record generation 0 as
                                      ancient samples. Must be `True` for
                                      tree sequence "recapitation". See
                                      :ref:`finishwithmsprime`.
    :type preserve_first_generation: bool
    :param check_demographic_event_timings: (True) If ``True``, then issue
                                            warnings if demographic events
                                            will occur prior to the current
                                            generation of the population.
    :type check_demographic_event_timings: bool

    The recording of genetic values into :attr:`fwdpy11.PopulationBase.genetic_values`
    is suppressed by default.  First, it is redundant with
    :attr:`fwdpy11.DiploidMetadata.g` for the common case of mutational effects on a
    single trait. Second, we save some memory by not tracking these matrices.
    However, it is useful to track these data for some
    cases when simulating multivariate mutational effects (pleiotropy).

    For a detailed description of `post_simplification_recorder`,
    see :ref:`tstimeseries`.

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

    """
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

    from ._fwdpy11 import MutationRegions
    from ._fwdpy11 import dispatch_create_GeneticMap
    from ._fwdpy11 import evolve_with_tree_sequences

    try:
        # DemographicModelDetails ?
        demographic_model = params.demography.model
    except AttributeError:
        # No, assume it is a DiscreteDemography
        demographic_model = params.demography
    except Exception as e:
        raise e

    if check_demographic_event_timings:
        _validate_event_timings(demographic_model, pop.generation)

    pneutral = 0.0
    if params.rates.neutral_mutation_rate + params.rates.selected_mutation_rate > 0.0:
        pneutral = params.rates.neutral_mutation_rate / (
            params.rates.neutral_mutation_rate + params.rates.selected_mutation_rate
        )
    mm = MutationRegions.create(pneutral, params.nregions, params.sregions)
    rm = dispatch_create_GeneticMap(params.rates.recombination_rate, params.recregions)

    from ._fwdpy11 import SampleRecorder

    sr = SampleRecorder()
    from ._fwdpy11 import _dgvalue_pointer_vector

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
        params.pself,
        params.prune_selected is False,
        suppress_table_indexing,
        record_gvalue_matrix,
        track_mutation_counts,
        remove_extinct_variants,
        reset_treeseqs_after_simplify,
        preserve_first_generation,
        post_simplification_recorder,
    )
