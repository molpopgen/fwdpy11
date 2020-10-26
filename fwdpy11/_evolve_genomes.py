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


def evolve_genomes(rng, pop, params, recorder=None):
    """
    Evolve a population without tree sequence recordings.  In other words,
    complete genomes must be simulated and tracked.

    :param rng: random number generator
    :type rng: :class:`fwdpy11.GSLrng`
    :param pop: A population
    :type pop: :class:`fwdpy11.DiploidPopulation`
    :param params: simulation parameters
    :type params: :class:`fwdpy11.ModelParams`
    :param recorder: (None) A temporal sampler/data recorder.
    :type recorder: callable

    .. note::
        If recorder is None,
        then :class:`fwdpy11.RecordNothing` will be used.

    .. versionchanged:: 0.8.0

        Update to refactored ModelParams
    """
    import warnings

    # warnings.simplefilter("default")
    warnings.warn(
        "Simulation without tree sequences is deprecated and will be removed in 0.11 or soon thereafter.", FutureWarning,
    )

    from ._fwdpy11 import MutationRegions
    from ._fwdpy11 import evolve_without_tree_sequences
    from ._fwdpy11 import dispatch_create_GeneticMap

    pneutral = 0.0
    if params.rates.neutral_mutation_rate + params.rates.selected_mutation_rate > 0.0:
        pneutral = params.rates.neutral_mutation_rate / (
            params.rates.neutral_mutation_rate + params.rates.selected_mutation_rate
        )
    mm = MutationRegions.create(pneutral, params.nregions, params.sregions)
    rm = dispatch_create_GeneticMap(params.rates.recombination_rate, params.recregions)

    if recorder is None:
        from ._fwdpy11 import RecordNothing

        recorder = RecordNothing()

    # FIXME: passing params.popsizes is a hack introduced in 0.6.0
    # top maintain API compatibility
    evolve_without_tree_sequences(
        rng,
        pop,
        params.popsizes,
        params.rates.neutral_mutation_rate,
        params.rates.selected_mutation_rate,
        mm,
        rm,
        params.gvalue,
        recorder,
        params.pself,
        params.prune_selected,
    )
