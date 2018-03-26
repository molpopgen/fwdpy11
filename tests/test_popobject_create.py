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


class testSlocusPopCreate(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.mutations = fwdpy11.VecMutation()
        self.fixations = fwdpy11.VecMutation()
        self.gametes = fwdpy11.VecGamete()
        self.diploids = fwdpy11.VecDiploid()
        self.mutations.append(fwdpy11.Mutation(0.1, -0.01, 1.0, 0, 0))
        self.fixations.append(fwdpy11.Mutation(0.1, -0.01, 1.0, 0, 0))
        self.gametes.append(fwdpy11.Gamete(
            (2, fwdpy11.VecUint32([]), fwdpy11.VecUint32([0]))))
        self.diploids.append(fwdpy11.SingleLocusDiploid(0, 0))

    def testConstruction(self):
        pop = fwdpy11.SlocusPop(self.diploids, self.gametes, self.mutations)
        self.assertTrue(pop.N, 1)

    def testStaticMethod(self):
        pop = fwdpy11.SlocusPop.create(
            self.diploids, self.gametes, self.mutations)
        self.assertTrue(type(pop) is fwdpy11.SlocusPop)
        # Test that data were moved and not copied:
        self.assertEqual(len(self.diploids), 0)
        self.assertEqual(len(self.gametes), 0)
        self.assertEqual(len(self.mutations), 0)

    def testStaticMethodWithFixations(self):
        ftimes = fwdpy11.VecUint32([1])
        pop = fwdpy11.SlocusPop.create(self.diploids,
                                       self.gametes,
                                       self.mutations,
                                       self.fixations,
                                       ftimes, 2)
        self.assertTrue(type(pop) is fwdpy11.SlocusPop)
        self.assertEqual(len(pop.fixations), len(pop.fixation_times))
        # Test that data were moved and not copied:
        self.assertEqual(len(self.diploids), 0)
        self.assertEqual(len(self.gametes), 0)
        self.assertEqual(len(self.mutations), 0)
        self.assertEqual(len(self.fixations), 0)
        self.assertEqual(len(ftimes), 0)


# class testSlocusPopGeneralMutVecCreate(unittest.TestCase):
#     @classmethod
#     def setUp(self):
#         self.mutations = fwdpy11.VecGeneralMutVec()
#         self.fixations = fwdpy11.VecGeneralMutVec()
#         self.gametes = fwdpy11.VecGamete()
#         self.diploids = fwdpy11.VecDiploid()
#         s = fwdpy11.VecDouble([-0.1, 0.1])
#         h = fwdpy11.VecDouble([0.0, 1.0])
#         m = fwdpy11.GeneralMutVec((s, h, 0.1, 3, 0))
#         self.mutations.append(m)
#         self.fixations.append(m)
#         self.gametes.append(fwdpy11.Gamete(
#             (2, fwdpy11.VecUint32([]), fwdpy11.VecUint32([0]))))
#         self.diploids.append(fwdpy11.SingleLocusDiploid(0, 0))
# 
#     def testConstruction(self):
#         pop = fwdpy11.SlocusPopGeneralMutVec(
#             self.diploids, self.gametes, self.mutations)
#         self.assertTrue(pop.N, 1)
# 
#     def testStaticMethod(self):
#         pop = fwdpy11.SlocusPopGeneralMutVec.create(
#             self.diploids, self.gametes, self.mutations)
#         self.assertTrue(type(pop) is fwdpy11.SlocusPopGeneralMutVec)
#         # Test that data were moved and not copied:
#         self.assertEqual(len(self.diploids), 0)
#         self.assertEqual(len(self.gametes), 0)
#         self.assertEqual(len(self.mutations), 0)
# 
#     def testStaticMethodWithFixations(self):
#         ftimes = fwdpy11.VecUint32([1])
#         pop = fwdpy11.SlocusPopGeneralMutVec.create(self.diploids,
#                                                     self.gametes,
#                                                     self.mutations,
#                                                     self.fixations,
#                                                     ftimes, 2)
#         self.assertTrue(type(pop) is fwdpy11.SlocusPopGeneralMutVec)
#         self.assertEqual(len(pop.fixations), len(pop.fixation_times))
#         # Test that data were moved and not copied:
#         self.assertEqual(len(self.diploids), 0)
#         self.assertEqual(len(self.gametes), 0)
#         self.assertEqual(len(self.mutations), 0)
#         self.assertEqual(len(self.fixations), 0)
#         self.assertEqual(len(ftimes), 0)


class testMlocusPopCreate(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.mutations = fwdpy11.VecMutation()
        self.fixations = fwdpy11.VecMutation()
        self.gametes = fwdpy11.VecGamete()
        self.diploids = fwdpy11.VecVecDiploid()
        self.mutations.append(fwdpy11.Mutation(0.1, -0.01, 1.0, 0, 0))
        self.fixations.append(fwdpy11.Mutation(0.1, -0.01, 1.0, 0, 0))
        self.gametes.append(fwdpy11.Gamete(
            (4, fwdpy11.VecUint32([]), fwdpy11.VecUint32([0]))))
        self.diploids.append(fwdpy11.VecDiploid(
            [fwdpy11.SingleLocusDiploid(0, 0)] * 2))

    def testConstruction(self):
        pop = fwdpy11.MlocusPop(self.diploids, self.gametes, self.mutations)
        self.assertTrue(pop.N, 1)

    def testStaticMethod(self):
        pop = fwdpy11.MlocusPop.create(
            self.diploids, self.gametes, self.mutations)
        self.assertTrue(type(pop) is fwdpy11.MlocusPop)
        # Test that data were moved and not copied:
        self.assertEqual(len(self.diploids), 0)
        self.assertEqual(len(self.gametes), 0)
        self.assertEqual(len(self.mutations), 0)

    def testStaticMethodWithFixations(self):
        ftimes = fwdpy11.VecUint32([1])
        pop = fwdpy11.MlocusPop.create(self.diploids,
                                       self.gametes,
                                       self.mutations,
                                       self.fixations,
                                       ftimes, 2)
        self.assertTrue(type(pop) is fwdpy11.MlocusPop)
        self.assertEqual(len(pop.fixations), len(pop.fixation_times))
        # Test that data were moved and not copied:
        self.assertEqual(len(self.diploids), 0)
        self.assertEqual(len(self.gametes), 0)
        self.assertEqual(len(self.mutations), 0)
        self.assertEqual(len(self.fixations), 0)
        self.assertEqual(len(ftimes), 0)


if __name__ == "__main__":
    unittest.main()
