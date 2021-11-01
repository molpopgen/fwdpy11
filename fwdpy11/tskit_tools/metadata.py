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

from .._fwdpy11 import Mutation
from ._flags import (INDIVIDUAL_IS_ALIVE, INDIVIDUAL_IS_FIRST_GENERATION,
                     INDIVIDUAL_IS_PRESERVED)


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
    alive: bool = attr.ib(converter=bool)
    preserved: bool = attr.ib(converter=bool)
    first_generation: bool = attr.ib(converter=bool)
    parents: typing.List[int]
    geography: typing.List[float]
    nodes: typing.List[int]

    @classmethod
    def from_table_row(cls, ind: tskit.IndividualTableRow):
        """
        Generate a metadata object from a :class:`tskit.IndividualTableRow`

        .. versionadded:: 0.15.0
        """
        alive = ind.flags & INDIVIDUAL_IS_ALIVE
        preserved = ind.flags & INDIVIDUAL_IS_PRESERVED
        first_generation = ind.flags & INDIVIDUAL_IS_FIRST_GENERATION
        return cls(
            **ind.metadata,
            alive=alive,
            preserved=preserved,
            first_generation=first_generation
        )


def decode_individual_metadata(
    ts: tskit.TreeSequence, rows: typing.Optional[typing.Union[int, slice]] = None
) -> typing.List[DiploidMetadata]:
    """
    Decodes a :class:`tskit.IndividualTable`.

    :param ts: A tree sequence
    :type ts: :class:`tskit.TreeSequence`

    :param rows: The rows to decode. If `None`, the entire table is decoded and returned.
    :type rows: int or slice or numpy array

    :returns: individual metadata
    :rtype: List of :class:`fwdpy11.tskit_tools.DiploidMetadata`

    .. warning::

        Decoding the entire table may consume a lot of system memory.

    .. versionadded:: 0.12.0

    .. versionchanged:: 0.15.0

        Add index/slice access to the table.

    .. versionchanged:: 0.17.0

        Change input from TableCollection to TreeSequence
    """
    rv = []

    if rows is None:
        _rows = slice(0, ts.num_individuals, 1)
    else:
        _rows = rows

    try:
        ind = ts.individual(_rows)
        rv.append(DiploidMetadata.from_table_row(ind))
    except Exception as _:
        try:
            for i in _rows:
                ind = ts.individual(i)
                rv.append(DiploidMetadata.from_table_row(ind))
        except Exception as _:
            for i in range(_rows.start, _rows.stop, _rows.step):
                ind = ts.individual(i)
                rv.append(DiploidMetadata.from_table_row(ind))

    return rv


def _append_mutation_metadata(ts, site, md, mutations):
    if md is not None:
        if "esizes" not in md:
            mutations.append(
                Mutation(
                    pos=ts.site(site).position,
                    s=md["s"],
                    h=md["h"],
                    g=md["origin"],
                    label=md["label"],
                )
            )
        else:
            mutations.append(
                Mutation(
                    pos=ts.site(site).position,
                    s=md["s"],
                    h=md["h"],
                    g=md["origin"],
                    esizes=md["esizes"],
                    heffects=md["heffects"],
                    label=md["label"],
                )
            )
    else:
        mutations.append(None)


def decode_mutation_metadata(
    ts: tskit.TreeSequence, rows: typing.Optional[typing.Union[int, slice]] = None
) -> typing.List[typing.Optional[Mutation]]:
    """
    Decodes metadata from a :class:`tskit.MutationTable`.

    :param ts: A tree sequence
    :type ts: :class:`tskit.TreeSequence`
    :param rows: The rows to decode. If `None`, the entire table is decoded and returned.
    :type rows: int or slice or numpy array

    :returns: Mutations
    :rtype: list

    .. warning::

        Decoding the entire table may consume a lot of system memory.

    .. versionadded:: 0.12.0

    .. versionchanged:: 0.15.0

        Return type is now a list, allowing
        some elements to be `None`.

        Add index/slice access to the table.

    .. versionchanged:: 0.17.0

        Change input from TableCollection to TreeSequence
    """
    mutations = []
    if rows is None:
        _rows = slice(0, ts.num_mutations, 1)
    else:
        _rows = rows

    try:
        m = ts.mutation(_rows)
        _append_mutation_metadata(ts, m.site, m.metadata, mutations)
    except Exception as _:
        try:
            for i in _rows:
                m = ts.mutation(i)
                _append_mutation_metadata(ts, m.site, m.metadata, mutations)
        except Exception as _:
            for i in range(_rows.start, _rows.stop, _rows.step):
                m = ts.mutation(i)
                _append_mutation_metadata(ts, m.site, m.metadata, mutations)
    return mutations
