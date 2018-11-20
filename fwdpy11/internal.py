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


def makeMutationRegions(rng, pop, neutral, selected, pneutral):
    """
    Convert user input into :class:`~fwdpy11.fwdpp_extensions.MutationRegions`

    :param rng: A :class:`fwdpy11.GSLrng`
    :param pop: A :class:`fwdp11.Population`
    :param neutral: A list of :class:`fwdpy11.regions.Region` objects.
    :param selected: A list of :class:`fwdpy11.regions.Sregion` objects.
    :param pneutral: The probability that a new mutation is neutral
    :ptype neutral: float

    :rtype: :class:`fwdpy11.fwdpp_extensions.MutationRegions`

    .. note:: Used by various "evolve" functions.
        Users probably won't need to call this.

        >>> import fwdpy11 as fp11
        >>> nregions = [fp11.Region(0,0.5,1),fp11.Region(1,1.5,1)]
        >>> sregions = [fp11.ExpS(0,0.5,1,1),fp11.GaussianS(1,1.5,1,0.25)]
        >>> mr = fp11.makeMutationRegions(nregions,sregions,0.1)
        >>> type(mr)
        <class 'fwdpy11.fwdpp_extensions.MutationRegions'>

    One or both lists may be empty:

        >>> mr = fp11.makeMutationRegions([],[],0.0)

    Neither list may be None:

        >>> mr = fp11.makeMutationRegions([],None,0.0)
        Traceback (most recent call last):
            ...
        TypeError: 'NoneType' object is not iterable

    .. versionchanged:: 0.2.0
        Added pneutral to handle changes in fwdpp 0.6.0
    """
    import math
    if math.isfinite(pneutral) is False:
        raise ValueError("pneutral not finite")
    elif pneutral < 0.0 or pneutral > 1.0:
        raise ValueError("pneutral must be in the range [0.0, 1.0]")
    import numpy as np

    # We need to reweight the user input so that new
    # mutations come out at the expected rates.
    # This is necessary b/c the user inputs weights
    # *separately* for neutral and selected variants.
    # These weights have to get combined into a single
    # weight vector on the C++ side, meaning that we
    # have some normalization to do:
    neutral_region_weights = np.array([i.w for i in neutral], dtype=np.float64)
    neutral_region_weights /= neutral_region_weights.sum()
    neutral_region_weights *= pneutral
    selected_region_weights = np.array([i.w for i in selected], dtype=np.float64)
    selected_region_weights /= selected_region_weights.sum()
    selected_region_weights *= (1.0-pneutral)
    summed_weights = neutral_region_weights.sum() + selected_region_weights.sum()
    neutral_region_weights /= summed_weights
    selected_region_weights /= summed_weights
    ntuples = [(i.b, i.e, j, i.l)
               for i, j in zip(neutral, neutral_region_weights)]
    stuples = [(i.b, i.e, j, i.l)
               for i, j in zip(selected, selected_region_weights)]
    callbacks = [i.callback() for i in selected]
    from .fwdpp_extensions import MutationRegions
    return MutationRegions(rng, pop, ntuples, stuples, callbacks)


def makeRecombinationRegions(rng, recrate, regions):
    """
    Convert user input into
    :class:`~fwdpy11.fwdpp_extensions.RecombinationRegions`

    :param rng: :class:`fwdpy11.GSLrng`
    :param recrate: (float) Recombination reate.
    :param regions: A list of :class:`fwdpy11.Region`

    :rtype: :class:`fwdpy11.fwdpp_extensions.RecombinationRegions`


    .. note::
        Used by various "evolve" functions.
        Users probably won't need to call this.

        >>> import fwdpy11 as fp11
        >>> rng = fwdpy11.GSLrng(42)
        >>> rregions = [fp11.Region(0,0.5,1),fp11.Region(1,1.5,1)]
        >>> rr = fp11.makeRecombinationRegions(rng, 1e-3, rregions)
        >>> type(rr)
        <class 'fwdpy11.fwdpp_extensions.RecombinationRegions'>

    """
    beg = [i.b for i in regions]
    end = [i.e for i in regions]
    weights = [i.w for i in regions]
    from .fwdpp_extensions import RecombinationRegions
    return RecombinationRegions(rng, recrate, beg, end, weights)
