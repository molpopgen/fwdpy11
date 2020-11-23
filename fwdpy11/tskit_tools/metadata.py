#
# Copyright (C) 2020 Kevin Thornton <krthornt@uci.edu>
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
import typing

import attr
import tskit

import fwdpy11


@attr.s(auto_attribs=True, frozen=True, kw_only=True, eq=True)
class DiploidMetadata(object):
    """
    Instances of this type are returned by
    :func:`fwdpy11.tskit_tools.decode_individual_metadata`.

    The following ``kwargs`` are used to construct instances of
    this type. Member attributes have the same names as these
    ``kwargs```:

    :param g: Genetic value
    :type g: float
    :param e: "Random/environmental" value
    :type e: float
    :param w: Fitness
    :type w: float
    :param sex: Sex
    :type sex: int
    :param deme: Deme
    :type deme: int
    :param label: Index in the population at time of living/sampling
    :type label: int
    :param alive: Individual is alive
    :type alive: bool
    :param preserved: Individual is preserved
    :type preserved: bool
    :param first_generation: Individual is first generation
    :type first_generation: bool
    :param parents: ``label`` of parents
    :type parents: list
    :param geography: Location in 3d space
    :type geography: list
    :param nodes: Nodes on tree sequence
    :type nodes: list

    .. versionadded:: 0.12.0
    """

    g: float
    e: float
    w: float
    sex: int
    deme: int
    label: int
    alive: bool
    preserved: bool
    first_generation: bool
    parents: typing.List[int]
    geography: typing.List[float]
    nodes: typing.List[int]


def decode_individual_metadata(tc: tskit.TableCollection):
    """
    Decodes a :class:`tskit.IndividualTable`.

    :param tc: A table collection
    :type tc: :class:`tskit.TableCollection`

    :returns: individual metadata
    :rtype: List of :class:`fwdpy11.tskit_tools.DiploidMetadata`


    .. versionadded:: 0.12.0
    """
    rv = []
    for i in tc.individuals:
        alive = i.flags & fwdpy11.tskit_tools.INDIVIDUAL_IS_ALIVE
        preserved = i.flags & fwdpy11.tskit_tools.INDIVIDUAL_IS_PRESERVED
        first_generation = i.flags & fwdpy11.tskit_tools.INDIVIDUAL_IS_FIRST_GENERATION
        rv.append(
            DiploidMetadata(
                **i.metadata,
                alive=alive,
                preserved=preserved,
                first_generation=first_generation
            )
        )

    return rv


def decode_mutation_metadata(tc: tskit.TableCollection) -> fwdpy11.MutationVector:
    """
    Decodes a :class:`tskit.MutationTable`.

    :param tc: A table collection
    :type tc: :class:`tskit.TableCollection`

    :returns: Mutations
    :rtype: :class:`fwdpy11.MutationVector`

    .. versionadded:: 0.12.0
    """
    mutations = fwdpy11.MutationVector()
    for m in tc.mutations:
        md = m.metadata
        if "esizes" not in md:
            mutations.append(
                fwdpy11.Mutation(
                    pos=tc.sites.position[m.site],
                    s=md["s"],
                    h=md["h"],
                    g=md["origin"],
                    label=md["label"],
                )
            )
        else:
            mutations.append(
                fwdpy11.Mutation(
                    pos=tc.sites.position[m.site],
                    s=md["s"],
                    h=md["h"],
                    g=md["origin"],
                    esizes=md["esizes"],
                    heffects=md["heffects"],
                    label=md["label"],
                )
            )
    return mutations
