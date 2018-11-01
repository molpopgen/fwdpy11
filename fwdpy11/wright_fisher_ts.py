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


def evolve(rng, pop, params, simplification_interval, recorder=None):
    """
    Evolve a population

    :param rng: An instance of :class:`fwdpy11.GSLrng`
    :param pop: An instance of :class:`fwdpy11.SlocusPop`
    :param params: An instance of :class:`fwdpy11.model_params.SlocusParams`
    :param simplification_interval: Number of generations between simplifications.
    :param recorder: (None) A temporal sampler/data recorder.

    .. note::
        If recorder is None,
        then :class:`fwdpy11.temporal_samplers.RecordNothing` will be used.

    """
    # import fwdpy11.SlocusPop
    # import fwdpy11.MlocusPop
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
    # if pop.__class__ is fwdpy11.SlocusPop:
    from .wright_fisher_slocus_ts import WFSlocusPop_ts
    # TODO: update to allow neutral mutations
    pneutral = 0  # params.mutrate_n/(params.mutrate_n+params.mutrate_s)
    mm = makeMutationRegions(rng, pop, params.nregions,
                             params.sregions, pneutral)
    rm = makeRecombinationRegions(rng, params.recrate, params.recregions)

    if recorder is None:
        from fwdpy11._tsevolveutils import NoAncientSamples
        recorder = NoAncientSamples()

    from fwdpy11._tsevolveutils import SampleRecorder
    sr = SampleRecorder()
    WFSlocusPop_ts(rng, pop, sr, simplification_interval, params.demography, params.mutrate_s,
                   params.recrate, mm, rm, params.gvalue,
                   recorder, params.pself, params.prune_selected is True)
    # else
    #     from .wright_fisher_mlocus import WFMlocusPop
    #     mm = [makeMutationRegions(rng, pop, i, j, n/(n+s)) for
    #           i, j, n, s in zip(params.nregions,
    #                             params.sregions,
    #                             params.mutrates_n, params.mutrates_s)]
    #     rm = [makeRecombinationRegions(rng, i, j) for i, j in zip(
    #         params.recrates, params.recregions)]

    #     if recorder is None:
    #         from fwdpy11.temporal_samplers import RecordNothing
    #         recorder = RecordNothing()

    #     WFMlocusPop(rng, pop, params.demography, params.mutrates_n,
    #                 params.mutrates_s, mm, rm, params.interlocus_rec,
    #                 params.gvalue, recorder, params.pself,
    #                 params.prune_selected)
