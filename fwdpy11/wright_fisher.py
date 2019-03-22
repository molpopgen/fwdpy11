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


def evolve(rng, pop, params, recorder=None):
    """
    Evolve a population

    :param rng: random number generator
    :type rng: :class:`fwdpy11.GSLrng`
    :param pop: A population
    :type pop: :class:`fwdpy11.DiploidPopulation`
    :param params: simulation parameters
    :type params: :class:`fwdpy11.model_params.ModelParams`
    :param recorder: (None) A temporal sampler/data recorder.
    :type recorder: callable

    .. note::
        If recorder is None,
        then :class:`fwdpy11.temporal_samplers.RecordNothing` will be used.

    """
    import fwdpy11.DiploidPopulation
    import warnings
    # Test parameters while suppressing warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # Will throw exception if anything is wrong:
        params.validate()

    from fwdpy11 import MutationRegions
    from fwdpy11 import RecombinationRegions
    from .wright_fisher_slocus import WFDiploidPopulation
    from fwdpy11 import GeneralizedGeneticMap
    pneutral = params.mutrate_n/(params.mutrate_n+params.mutrate_s)
    mm = MutationRegions.create(pneutral, params.nregions, params.sregions)
    if all([i.__class__ is fwdpy11.Region for i in params.recregions]) is True:
        rm = RecombinationRegions(params.recrate, params.recregions)
    else:
        rm = GeneralizedGeneticMap(params.recregions)

    if recorder is None:
        from fwdpy11.temporal_samplers import RecordNothing
        recorder = RecordNothing()

    WFDiploidPopulation(rng, pop, params.demography,
                params.mutrate_n, params.mutrate_s,
                params.recrate, mm, rm, params.gvalue,
                recorder, params.pself, params.prune_selected)
