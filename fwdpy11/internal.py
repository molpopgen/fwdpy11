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
def makeMutationRegions(neutral,selected):
    """
    Convert user input into :class:`~fwdpy11.fwdpp_extensions.MutationRegions`

    :param neutral: A list of :class:`fwdpy11.regions.Region` objects.
    :param selected: A list of :class:`fwdpy11.regions.Sregion` objects.

    :rtype: :class:`fwdpy11.fwdpp_extensions.MutationRegions`

    .. note:: Used by various "evolve" functions.  Users probably won't need to call this.
    
        >>> import fwdpy11 as fp11
        >>> nregions = [fp11.Region(0,0.5,1),fp11.Region(1,1.5,1)]
        >>> sregions = [fp11.ExpS(0,0.5,1,1),fp11.GaussianS(1,1.5,1,0.25)]
        >>> mr = fp11.makeMutationRegions(nregions,sregions)
        >>> type(mr)
        <class 'fwdpy11.fwdpp_extensions.MutationRegions'>

    One or both lists may be empty:

        >>> mr = fp11.makeMutationRegions([],[])

    Neither list may be None:

        >>> mr = fp11.makeMutationRegions([],None)
        Traceback (most recent call last):
            ...
        TypeError: 'NoneType' object is not iterable
    """
    nbeg = [i.b for i in neutral]
    nend = [i.e for i in neutral]
    nweights = [i.w for i in neutral]
    sbeg = [i.b for i in selected]
    send = [i.e for i in selected]
    sweights = [i.w for i in selected]
    sh = [i.callback() for i in selected]
    from .fwdpp_extensions import MutationRegions
    return MutationRegions(nbeg,nend,nweights,sbeg,send,sweights,sh)

def makeRecombinationRegions(regions):
    """
    Convert user input into :class:`~fwdpy11.fwdpp_extensions.RecombinationRegions`

    :param neutral: A list of :class:`fwdpy11.regions.Region` objects.
    :param selected: A list of :class:`fwdpy11.regions.Sregion` objects.

    :rtype: :class:`fwdpy11.fwdpp_extensions.MutationRegions`

    .. note:: Used by various "evolve" functions.  Users probably won't need to call this.
    
        >>> import fwdpy11 as fp11
        >>> rregions = [fp11.Region(0,0.5,1),fp11.Region(1,1.5,1)]
        >>> rr = fp11.makeRecombinationRegions(rregions)
        >>> type(rr)
        <class 'fwdpy11.fwdpp_extensions.RecombinationRegions'>

    """
    beg = [i.b for i in regions]
    end = [i.e for i in regions]
    weights = [i.w for i in regions]
    from .fwdpp_extensions import RecombinationRegions
    return RecombinationRegions(beg,end,weights)
