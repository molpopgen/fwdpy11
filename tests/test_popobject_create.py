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

# TYPETEST

import fwdpy11
import unittest
import pickle


class testDiploidCreate(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.d = fwdpy11.DiploidGenotype(0, 0)

    def testPickle(self):
        x = pickle.dumps(self.d)
        dd = pickle.loads(x)
        self.assertEqual(dd.first, 0)
        self.assertEqual(dd.second, 0)


class testDiploidPopulationCreate(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.mutations = fwdpy11.MutationVector()
        self.fixations = fwdpy11.MutationVector()
        self.gametes = fwdpy11.GameteVector()
        self.diploids = fwdpy11.DiploidVector()
        self.mutations.append(fwdpy11.Mutation(0.1, -0.01, 1.0, 0, 0))
        self.fixations.append(fwdpy11.Mutation(0.1, -0.01, 1.0, 0, 0))
        self.gametes.append(fwdpy11.Gamete(
            (2, fwdpy11.VecUint32([]), fwdpy11.VecUint32([0]))))
        self.diploids.append(fwdpy11.DiploidGenotype(0, 0))

    def testConstruction(self):
        pop = fwdpy11.DiploidPopulation(self.diploids, self.gametes, self.mutations)
        self.assertTrue(pop.N, 1)

    def testStaticMethod(self):
        pop = fwdpy11.DiploidPopulation.create(
            self.diploids, self.gametes, self.mutations)
        self.assertTrue(type(pop) is fwdpy11.DiploidPopulation)
        # Test that data were moved and not copied:
        self.assertEqual(len(self.diploids), 0)
        self.assertEqual(len(self.gametes), 0)
        self.assertEqual(len(self.mutations), 0)

    def testStaticMethodWithFixations(self):
        ftimes = fwdpy11.VecUint32([1])
        pop = fwdpy11.DiploidPopulation.create(self.diploids,
                                       self.gametes,
                                       self.mutations,
                                       self.fixations,
                                       ftimes, 2)
        self.assertTrue(type(pop) is fwdpy11.DiploidPopulation)
        self.assertEqual(len(pop.fixations), len(pop.fixation_times))
        # Test that data were moved and not copied:
        self.assertEqual(len(self.diploids), 0)
        self.assertEqual(len(self.gametes), 0)
        self.assertEqual(len(self.mutations), 0)
        self.assertEqual(len(self.fixations), 0)
        self.assertEqual(len(ftimes), 0)


class testGameteKeySorting(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.mutations = fwdpy11.MutationVector()
        self.fixations = fwdpy11.MutationVector()
        self.gametes = fwdpy11.GameteVector()
        self.diploids = fwdpy11.DiploidVector()
        self.mutations.append(fwdpy11.Mutation(0.1, -0.01, 1.0, 0, 0))
        self.mutations.append(fwdpy11.Mutation(0.01, -2, 1, 0, 0))
        self.gametes.append(fwdpy11.Gamete(
            (1, fwdpy11.VecUint32([]), fwdpy11.VecUint32([0]))))
        self.gametes.append(fwdpy11.Gamete(
            (1, fwdpy11.VecUint32([]), fwdpy11.VecUint32([0, 1]))))
        self.diploids.append(fwdpy11.DiploidGenotype(0, 1))

    def testUnsortedGametes(self):
        with self.assertRaises(ValueError):
            pop = fwdpy11.DiploidPopulation.create(
                self.diploids, self.gametes, self.mutations)

    def testSortingGametes(self):
        fwdpy11.sort_gamete_keys(self.gametes, self.mutations)
        pop = fwdpy11.DiploidPopulation.create(
            self.diploids, self.gametes, self.mutations)


if __name__ == "__main__":
    unittest.main()
