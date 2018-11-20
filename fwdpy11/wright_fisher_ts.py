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
           suppress_table_indexing=False):
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

    .. note::
        If recorder is None,
        then :class:`fwdpy11.tsrecorders.NoAncientSamples` will be used.

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

    from .internal import makeMutationRegions, makeRecombinationRegions
    from .wright_fisher_slocus_ts import WFSlocusPop_ts
    # TODO: update to allow neutral mutations
    pneutral = 0
    mm = makeMutationRegions(rng, pop, params.nregions,
                             params.sregions, pneutral)
    rm = makeRecombinationRegions(rng, params.recrate, params.recregions)

    if recorder is None:
        from fwdpy11.tsrecorders import NoAncientSamples
        recorder = NoAncientSamples()

    from fwdpy11.tsrecorders import SampleRecorder
    sr = SampleRecorder()
    WFSlocusPop_ts(rng, pop, sr, simplification_interval,
                   params.demography, params.mutrate_s,
                   params.recrate, mm, rm, params.gvalue,
                   recorder, params.pself, params.prune_selected is True,
                   suppress_table_indexing)
