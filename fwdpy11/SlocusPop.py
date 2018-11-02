#
# Copyright(C) 2017 Kevin Thornton < krthornt @uci.edu>
#
# This file is part of fwdpy11.
#
# fwdpy11 is free software : you can redistribute it and / or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# fwdpy11 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with fwdpy11.If not, see < http: //www.gnu.org/licenses/>.
#
from fwdpy11._Populations import _SlocusPop


class SlocusPop(_SlocusPop):
    """
    Population object representing a single
    locus.
    """
    pass

    @staticmethod
    def create(diploids, gametes, mutations, *args):
        """
        Create a new object from input data.
        Unlike the constructor method, this method results
        in no temporary copies of input data.

        :param diplods: A :class:`fwdpy11.VecDiploid`
        :param gametes: A :class:`fwdpy11.VecGamete`
        :param mutations: A :class:`fwdpy11.VecMutation`
        :param args: Fixations, fixation times, and generation

        :rtype: :class:`fwdpy11.SlocusPop`

        See :ref:`popobjects` for example use.

        When passing in extra args, they must be the following:

        fixations: A :class:`fwdpy11.VecMutation`
        fixation times: A :class:`fwdpy11.VecUint32`
        generation: A non-negative integer

        It is required that len(fixations) == len(fixation times).

        The result of passing in these extra args will be an object
        with its fixation data populated and its generation set
        to the input value.

        .. versionadded:: 0.1.4

        """
        return SlocusPop(_SlocusPop.create(diploids,
                                           gametes, mutations,
                                           args))

    @staticmethod
    def create_from_msprime(ts, discretize_time):
        from fwdpy11.ts_from_msprime import _create_SlocusPop
        return _create_SlocusPop(ts, discretize_time)
