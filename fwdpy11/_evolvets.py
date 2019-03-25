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


def evolvets(rng, pop, params, simplification_interval, recorder=None,
           suppress_table_indexing=False, record_gvalue_matrix=False,
           stopping_criterion=None,
           track_mutation_counts=False,
           remove_extinct_variants=True):
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
        then :class:`fwdpy11.NoAncientSamples` will be used.

    """
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
        from ._fwdpy11 import NoAncientSamples
        recorder = NoAncientSamples()

    if stopping_criterion is None:
        from ._fwdpy11 import _no_stopping
        stopping_criterion = _no_stopping

    from ._fwdpy11 import MutationRegions
    from ._fwdpy11 import dispatch_create_GeneticMap
    from ._fwdpy11 import evolve_with_tree_sequences
    # TODO: update to allow neutral mutations
    pneutral = 0
    mm = MutationRegions.create(pneutral, params.nregions, params.sregions)
    rm = dispatch_create_GeneticMap(params.recrate, params.recregions)

    from ._fwdpy11 import SampleRecorder
    sr = SampleRecorder()
    evolve_with_tree_sequences(rng, pop, sr, simplification_interval,
                               params.demography, params.mutrate_s,
                               mm, rm, params.gvalue,
                               recorder, stopping_criterion,
                               params.pself, params.prune_selected is False,
                               suppress_table_indexing, record_gvalue_matrix,
                               track_mutation_counts,
                               remove_extinct_variants)
