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


def evolve(rng, pop, params, simplification_interval, recorder=None,
           suppress_table_indexing=False, record_gvalue_matrix=False,
           stopping_criterion=None,
           track_mutation_counts=False,
           remove_extinct_variants=True):
    """
    Evolve a population

    :param rng: random number generator
    :type rng: :class:`fwdpy11.GSLrng`
    :param pop: A population
    :type pop: :class:`fwdpy11.SlocusPop`
    :param params: simulation parameters
    :type params: :class:`fwdpy11.model_params.ModelParams`
    :param simplification_interval: Number of generations between simplifications.
    :type simplification_interval: int
    :param recorder: (None) A temporal sampler/data recorder.
    :type recorder: callable
    :param suppress_table_indexing: (False) Prevents edge table indexing until end of simulation
    :type suppress_table_indexing: boolean
    :param record_gvalue_matrix: (False) Whether to record genetic values into :attr:`fwdpy11.Population.genetic_values`.
    :type record_gvalue_matrix: boolean

    The recording of genetic values into :attr:`fwdpy11.Population.genetic_values` is supprssed by default.  First, it
    is redundant with :attr:`fwdpy11.DiploidMetadata.g` for the common case of mutational effects on a single trait.
    Second, we save some memory by not tracking these matrices.  However, it is useful to track these data for some
    cases when simulating multivariate mutational effects (pleiotropy).

    .. note::
        If recorder is None,
        then :class:`fwdpy11.tsrecorders.NoAncientSamples` will be used.

    """
    import fwdpy11.SlocusPop
    import fwdpy11.MlocusPop
    import warnings

    # Currently, we do not support simulating neutral mutations
    # during tree sequence simulations, so we make sure that there
    # are no neutral regions/rates:
    if len(params.nregions) != 0:
        raise ValueError(
            "Simulation of neutral mutations on tree sequences not supported (yet).")

    # Test parameters while suppressing warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # Will throw exception if anything is wrong:
        params.validate()

    if recorder is None:
        from fwdpy11.tsrecorders import NoAncientSamples
        recorder = NoAncientSamples()

    if stopping_criterion is None:
        from ._tsevolution import _no_stopping
        stopping_criterion = _no_stopping

    from fwdpy11 import MutationRegions
    from fwdpy11 import RecombinationRegions
    if pop.__class__ is fwdpy11.SlocusPop:
        from fwdpy11 import GeneralizedGeneticMap
        from ._tsevolution import WFSlocusPop_ts
        # TODO: update to allow neutral mutations
        pneutral = 0
        mm = MutationRegions.create(pneutral, params.nregions, params.sregions)
        if all([i.__class__ is fwdpy11.Region for i in params.recregions]) is True:
            rm = RecombinationRegions(params.recrate, params.recregions)
        else:
            rm = GeneralizedGeneticMap(params.recregions)

        from fwdpy11.tsrecorders import SampleRecorder
        sr = SampleRecorder()
        WFSlocusPop_ts(rng, pop, sr, simplification_interval,
                       params.demography, params.mutrate_s,
                       mm, rm, params.gvalue,
                       recorder, stopping_criterion,
                       params.pself, params.prune_selected is False,
                       suppress_table_indexing, record_gvalue_matrix,
                       track_mutation_counts,
                       remove_extinct_variants)
    else:
        from ._tsevolution import WFMlocusPop_ts
        from fwdpy11 import MlocusMutationRegions
        from fwdpy11 import MlocusRecombinationRegions
        mm = MlocusMutationRegions()
        # The user is allowed to input an empty list
        # for nregions.  However, to make our lives
        # sane, we need to detect that case
        # and create a proxy object of a list of empty
        # lists
        nregion_proxy = params.nregions
        if len(nregion_proxy) == 0:
            nregion_proxy = [[]]*len(params.sregions)
        for n, s, i, j in zip(params.mutrate_n,
                              params.mutrate_s,
                              nregion_proxy,
                              params.sregions):
            pn = n/(n+s)
            if pn != 0.0:
                raise ValueError("neutral mutations not allowed")
            temp = MutationRegions.create(pn, i, j)
            mm.append(temp)
        rm = MlocusRecombinationRegions()
        for i, j in zip(params.recrates, params.recregions):
            rm.append(RecombinationRegions(i, j))
        if recorder is None:
            from fwdpy11.temporal_samplers import RecordNothing
            recorder = RecordNothing()
        from fwdpy11.tsrecorders import SampleRecorder
        sr = SampleRecorder()
        WFMlocusPop_ts(rng, pop, sr, simplification_interval,
                       params.demography, params.mutrate_s,
                       mm, rm,
                       params.interlocus_rec, params.gvalue,
                       recorder, stopping_criterion,
                       params.pself, params.prune_selected is False,
                       suppress_table_indexing, record_gvalue_matrix,
                       track_mutation_counts,
                       remove_extinct_variants)
