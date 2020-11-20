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


@attr.s(auto_attribs=True, kw_only=True, eq=True)
class DiploidMetadata(object):
    g: float
    e: float
    w: float
    sex: int
    deme: int
    label: int
    parents: typing.List[int]
    geography: typing.List[float]
    nodes: typing.List[int]


def decode_individual_metadata(tc):
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
