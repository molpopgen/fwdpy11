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

    """
    import warnings
    # Test parameters while suppressing warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # Will throw exception if anything is wrong:
        params.validate()

    from ._fwdpy11 import MutationRegions
    from ._fwdpy11 import evolve_without_tree_sequences
    from ._fwdpy11 import dispatch_create_GeneticMap
    pneutral = params.mutrate_n/(params.mutrate_n+params.mutrate_s)
    mm = MutationRegions.create(pneutral, params.nregions, params.sregions)
    rm = dispatch_create_GeneticMap(params.recrate, params.recregions)

    if recorder is None:
        from ._fwdpy11 import RecordNothing
        recorder = RecordNothing()

    evolve_without_tree_sequences(rng, pop, params.demography,
                                  params.mutrate_n, params.mutrate_s,
                                  params.recrate, mm, rm, params.gvalue,
                                  recorder, params.pself, params.prune_selected)
