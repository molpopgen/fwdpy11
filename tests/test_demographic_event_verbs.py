#
# Copyright (C) 2017-2020 Kevin Thornton <krthornt@uci.edu>
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

"""
The concepts tested here are redundant with other test files.
The goal is to demonstrate specific examples that are hopefully
digestible and can serve as examples.
"""

import typing
import unittest

import attr
import numpy as np

import fwdpy11


@attr.s(auto_attribs=True)
class DemeSizeAtTime(object):
    generation: int
    deme: int
    size: int


@attr.s(auto_attribs=True)
class DemeSizeTracker(object):
    data: typing.List[DemeSizeAtTime]

    def __call__(
        self, pop: fwdpy11.DiploidPopulation, sampler: fwdpy11.SampleRecorder
    ) -> None:
        deme_sizes = pop.deme_sizes()
        for i, j in zip(deme_sizes[0], deme_sizes[1]):
            self.data.append(DemeSizeAtTime(pop.generation, i, j))


def branch(when: int, fraction: float) -> fwdpy11.DiscreteDemography:
    """
    A deme gives rise to a new daugther deme.
    The original deme continues to exist.
    """
    mass_migrations = [
        fwdpy11.copy_individuals(when=when, source=0, destination=1, fraction=fraction)
    ]
    return fwdpy11.DiscreteDemography(mass_migrations=mass_migrations)


def split(when: int, fraction: float) -> fwdpy11.DiscreteDemography:
    """
    Parental deme 0 gives rise to daughter demes 1 and 2.
    Deme 0 ceases to exist.

    Daughter demes are created by copying parents.
    """
    mass_migrations = [
        fwdpy11.copy_individuals(when=when, source=0, destination=1, fraction=fraction),
        fwdpy11.copy_individuals(
            when=when, source=0, destination=2, fraction=1.0 - fraction
        ),
    ]
    set_deme_sizes = [fwdpy11.SetDemeSize(when=when, deme=0, new_size=0)]
    return fwdpy11.DiscreteDemography(
        mass_migrations=mass_migrations, set_deme_sizes=set_deme_sizes
    )


def split_via_moves(when: int, fraction: float) -> fwdpy11.DiscreteDemography:
    """
    Parental deme 0 gives rise to daughter demes 1 and 2.
    Deme 0 ceases to exist.

    Daughter demes are created by moving parents.
    """
    mass_migrations = [
        fwdpy11.move_individuals(when=when, source=0, destination=1, fraction=fraction),
        fwdpy11.move_individuals(
            when=when, source=0, destination=2, fraction=1.0 - fraction
        ),
    ]
    return fwdpy11.DiscreteDemography(mass_migrations=mass_migrations)


def merge(
    when: int, p0: float, p1: float, keep_ancestor_demes: bool
) -> fwdpy11.DiscreteDemography:
    """
    Demes 0 and 1 merge to form deme 2, which is formed by fractions
    p0 and p1 of demes 0 and 1, respectively.

    Whether or not to "kill off" the ancestral demes is an input choice
    """
    mass_migrations = [
        fwdpy11.copy_individuals(when=when, source=0, destination=2, fraction=p0),
        fwdpy11.copy_individuals(when=when, source=1, destination=2, fraction=p1),
    ]
    if keep_ancestor_demes is False:
        set_deme_sizes = [
            fwdpy11.SetDemeSize(when=when, deme=0, new_size=0),
            fwdpy11.SetDemeSize(when=when, deme=1, new_size=0),
        ]
    else:
        set_deme_sizes = []
    return fwdpy11.DiscreteDemography(
        mass_migrations=mass_migrations, set_deme_sizes=set_deme_sizes
    )


def pulse_migration(when: int, fraction: float) -> fwdpy11.DiscreteDemography:
    """
    A fraction of individuals from deme 0 are copied to deme 1.

    If deme 1 does not exist at time when, it will be created.  Thus, this is the
    same as "branch".
    """
    return branch(when, fraction)


def build_params(
    demography: fwdpy11.DiscreteDemography, simlen: int
) -> fwdpy11.ModelParams:
    pdict = {
        "nregions": [],
        "sregions": [],
        "recregions": [],
        "rates": [0.0, 0.0, 0.0],
        "gvalue": fwdpy11.Multiplicative(2.0),
        "demography": demography,
        "simlen": simlen,
    }
    return fwdpy11.ModelParams(**pdict)


def run_model(
    params: fwdpy11.ModelParams, popsizes: typing.Union[int, typing.List], seed: int
) -> fwdpy11.DiploidPopulation:
    # 1.0 is the genome length for the tree sequence recording
    pop = fwdpy11.DiploidPopulation(popsizes, 1.0)
    rng = fwdpy11.GSLrng(seed)
    tracker = DemeSizeTracker([])
    fwdpy11.evolvets(rng, pop, params, 10, tracker)  # simplify every 10 generations
    return pop, tracker.data


class TestBranch(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.when = 10
        self.fraction = 0.25
        self.popsizes = 100
        demography = branch(self.when, self.fraction)
        self.params = build_params(demography, 50)
        self.pop, self.size_data = run_model(self.params, self.popsizes, 666)

    def test_final_sizes(self):
        deme_sizes = self.pop.deme_sizes(as_dict=True)
        self.assertEqual(deme_sizes[0], self.popsizes)
        self.assertEqual(deme_sizes[1], np.rint(self.fraction * self.popsizes))

    def test_sizes_over_time(self):
        self.assertEqual(
            len([i for i in self.size_data if i.deme == 0]), self.params.simlen
        )
        # Deme 1 first appeared in generation 10, so there were 9 generations
        # w/o that deme around.
        self.assertEqual(
            len([i for i in self.size_data if i.deme == 1]),
            self.params.simlen - self.when + 1,
        )


class TestSplit(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.when = 10
        self.fraction = 0.25
        self.popsizes = 100
        demography = split(self.when, self.fraction)
        self.params = build_params(demography, 50)
        self.pop, self.size_data = run_model(self.params, self.popsizes, 666)

    def test_final_sizes(self):
        deme_sizes = self.pop.deme_sizes(as_dict=True)
        self.assertTrue(0 not in deme_sizes)
        self.assertEqual(deme_sizes[1], np.rint(self.fraction * self.popsizes))
        self.assertEqual(deme_sizes[2], np.rint((1.0 - self.fraction) * self.popsizes))

    def test_sizes_over_time(self):
        self.assertEqual(
            max([i.generation for i in self.size_data if i.deme == 0]), self.when
        )
        self.assertEqual(
            min([i.generation for i in self.size_data if i.deme == 1]), self.when
        )
        self.assertEqual(
            min([i.generation for i in self.size_data if i.deme == 2]), self.when
        )


class TestSplitViaMoves(unittest.TestCase):
    # This model differs from the above in that
    # no parents will exist in deme 0 at time self.when - 1
    # because they were all moved out to create the other demes.
    @classmethod
    def setUpClass(self):
        self.when = 10
        self.fraction = 0.25
        self.popsizes = 100
        demography = split_via_moves(self.when, self.fraction)
        self.params = build_params(demography, 50)
        self.pop, self.size_data = run_model(self.params, self.popsizes, 666)

    def test_final_sizes(self):
        deme_sizes = self.pop.deme_sizes(as_dict=True)
        self.assertTrue(0 not in deme_sizes)
        self.assertEqual(deme_sizes[1], np.rint(self.fraction * self.popsizes))
        self.assertEqual(deme_sizes[2], np.rint((1.0 - self.fraction) * self.popsizes))

    def test_sizes_over_time(self):
        self.assertEqual(
            max([i.generation for i in self.size_data if i.deme == 0]), self.when - 1
        )
        self.assertEqual(
            min([i.generation for i in self.size_data if i.deme == 1]), self.when
        )
        self.assertEqual(
            min([i.generation for i in self.size_data if i.deme == 2]), self.when
        )


class TestMergeKeepAncestralDemes(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.when = 10
        self.p0 = 0.33
        self.p1 = 0.77
        self.popsizes = [100, 100]
        demography = merge(self.when, self.p0, self.p1, True)
        self.params = build_params(demography, 50)
        self.pop, self.size_data = run_model(self.params, self.popsizes, 666)

    def test_final_sizes(self):
        deme_sizes = self.pop.deme_sizes(as_dict=True)
        self.assertEqual(deme_sizes[0], self.popsizes[0])
        self.assertEqual(deme_sizes[1], self.popsizes[1])
        self.assertEqual(
            deme_sizes[2],
            np.rint(self.p0 * self.popsizes[0]) + np.rint(self.p1 * self.popsizes[1]),
        )

    def test_sizes_over_time(self):
        self.assertEqual(
            len([i for i in self.size_data if i.deme == 0]), self.params.simlen
        )
        self.assertEqual(
            len([i for i in self.size_data if i.deme == 1]), self.params.simlen
        )

        min_deme_2 = min([i.generation for i in self.size_data if i.deme == 2])

        self.assertEqual(min_deme_2, self.when)

        # Deme 2 first appeared in generation 10, so there were 9 generations
        # w/o that deme around.
        self.assertEqual(
            len([i for i in self.size_data if i.deme == 2]),
            self.params.simlen - self.when + 1,
        )


class TestMergeKillAncestralDemes(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.when = 10
        self.p0 = 0.33
        self.p1 = 0.77
        self.popsizes = [100, 100]
        demography = merge(self.when, self.p0, self.p1, False)
        self.params = build_params(demography, 50)
        self.pop, self.size_data = run_model(self.params, self.popsizes, 666)

    def test_final_sizes(self):
        deme_sizes = self.pop.deme_sizes(as_dict=True)
        self.assertTrue(0 not in deme_sizes)
        self.assertTrue(1 not in deme_sizes)
        self.assertEqual(
            deme_sizes[2],
            np.rint(self.p0 * self.popsizes[0]) + np.rint(self.p1 * self.popsizes[1]),
        )

    def test_sizes_over_time(self):
        self.assertEqual(len([i for i in self.size_data if i.deme == 0]), self.when)
        self.assertEqual(len([i for i in self.size_data if i.deme == 1]), self.when)

        min_deme_2 = min([i.generation for i in self.size_data if i.deme == 2])

        self.assertEqual(min_deme_2, self.when)

        # Deme 2 first appeared in generation 10, so there were 9 generations
        # w/o that deme around.
        self.assertEqual(
            len([i for i in self.size_data if i.deme == 2]),
            self.params.simlen - self.when + 1,
        )


class TestPulseMigration(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.when = 10
        self.fraction = 0.23
        self.popsizes = [200, 100]
        demography = pulse_migration(self.when, self.fraction)
        self.params = build_params(demography, 50)
        self.pop, self.size_data = run_model(self.params, self.popsizes, 666)

    def test_final_sizes(self):
        deme_sizes = self.pop.deme_sizes(as_dict=True)
        self.assertEqual(deme_sizes[0], self.popsizes[0])
        self.assertEqual(
            deme_sizes[1], np.rint(self.fraction * self.popsizes[0]) + self.popsizes[1]
        )

    def test_sizes_over_time(self):
        self.assertTrue(
            all([i.size == self.popsizes[0] for i in self.size_data if i.deme == 0])
        )
        self.assertTrue(
            all(
                [
                    i.size == self.popsizes[1]
                    for i in self.size_data
                    if i.deme == 1 and i.generation < 10
                ]
            )
        )
        self.assertTrue(
            all(
                [
                    i.size
                    == self.popsizes[1] + np.rint(self.fraction * self.popsizes[0])
                    for i in self.size_data
                    if i.deme == 1 and i.generation >= 10
                ]
            )
        )


if __name__ == "__main__":
    unittest.main()
