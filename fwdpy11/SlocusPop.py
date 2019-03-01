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
    def create_from_tskit(ts):
        """
        Create a new object from an tskit.TreeSequence

        :param ts: A tree sequence from tskit
        :type ts: tskit.TreeSequence

        :rtype: :class:`fwdpy11.SlocusPop`
        :returns: A population object with an initialized
        :class:`fwdpy11.ts.TableCollection`

        .. versionadded:: 0.2.0

        .. note::

            In general, initializing a population using
            the output from a coalescent simulation is
            a tricky business.  There are issues of
            parameter scaling and the appropriateness
            of the coalescent model itself. A key issue
            is that your input tree sequence must have
            node times in the correct time units! (Generations,
            for example.) See :ref:`ts` for more discussion

        """
        from fwdpy11.ts_from_tskit import _create_SlocusPop
        return SlocusPop(_create_SlocusPop(ts))

    def dump_tables_to_tskit(self):
        """
        Dump the population's TableCollection into
        an tskit TreeSequence

        :rtype: tskit.TreeSequence

        .. todo::

            Incorporate the various metadata values.
        """
        import fwdpy11._tables_to_tskit as t2msp
        return t2msp.dump_tables_to_tskit(self)
