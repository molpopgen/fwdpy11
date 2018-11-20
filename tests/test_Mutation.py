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

import unittest
import fwdpy11
import pickle


class testSingleEffectMutation(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.m = fwdpy11.Mutation((1.0, -1.0, 0.25, 0, 13))

    def testConstruct(self):
        self.assertEqual(self.m.pos, 1.0)
        self.assertEqual(self.m.s, -1.0)
        self.assertEqual(self.m.h, 0.25)
        self.assertEqual(self.m.g, 0)
        self.assertEqual(self.m.label, 13)
        self.assertEqual(len(self.m.esizes), 0)
        self.assertEqual(len(self.m.heffects), 0)

    def testPickle(self):
        p = pickle.dumps(self.m)
        up = pickle.loads(p)
        self.assertEqual(up, self.m)


class testMultiEffectMutation(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.m = fwdpy11.Mutation(
            1.0, -1.0, 0.25, 0, [1., 2.], [-1.0, 1. / 3.], 13)

    def testConstruct(self):
        self.assertEqual(self.m.pos, 1.0)
        self.assertEqual(self.m.s, -1.0)
        self.assertEqual(self.m.h, 0.25)
        self.assertEqual(self.m.g, 0)
        self.assertEqual(self.m.label, 13)
        self.assertEqual(len(self.m.esizes), 2)
        self.assertEqual(len(self.m.heffects), 2)

    def testPickle(self):
        p = pickle.dumps(self.m)
        up = pickle.loads(p)
        self.assertEqual(up, self.m)


class testMultiEffectMutationSlocusPopPickling(unittest.TestCase):
    @classmethod
    def setUp(self):
        mutations = fwdpy11.VecMutation()
        fixations = fwdpy11.VecMutation()
        gametes = fwdpy11.VecGamete()
        diploids = fwdpy11.VecDiploid()
        mutations.append(fwdpy11.Mutation(
            0.1, -0.01, 1.0, 0, [-1., 2.], [5., 4.], 0))
        fixations.append(fwdpy11.Mutation(
            0.1, -0.01, 1.0, 0, [-1., 2.], [5., 4.], 0))
        gametes.append(fwdpy11.Gamete(
            (2, fwdpy11.VecUint32([]), fwdpy11.VecUint32([0]))))
        diploids.append(fwdpy11.DiploidGenotype(0, 0))
        ftimes = fwdpy11.VecUint32([1])
        self.pop = fwdpy11.SlocusPop.create(diploids,
                                            gametes,
                                            mutations,
                                            fixations,
                                            ftimes, 2)

    def testPickle(self):
        p = pickle.dumps(self.pop)
        up = pickle.loads(p)
        self.assertEqual(self.pop, up)


class testMultiEffectMutationMlocusPopPickling(unittest.TestCase):
    @classmethod
    def setUp(self):
        mutations = fwdpy11.VecMutation()
        fixations = fwdpy11.VecMutation()
        gametes = fwdpy11.VecGamete()
        diploids = fwdpy11.VecVecDiploid()
        mutations.append(fwdpy11.Mutation(
            0.1, -0.01, 1.0, 0, [-1., 2.], [5., 4.], 0))
        fixations.append(fwdpy11.Mutation(
            0.1, -0.01, 1.0, 0, [-1., 2.], [5., 4.], 0))
        gametes.append(fwdpy11.Gamete(
            (4, fwdpy11.VecUint32([]), fwdpy11.VecUint32([0]))))
        diploids.append(fwdpy11.VecDiploid(
            [fwdpy11.DiploidGenotype(0, 0)] * 2))
        ftimes = fwdpy11.VecUint32([1])
        self.pop = fwdpy11.MlocusPop.create(diploids,
                                            gametes,
                                            mutations,
                                            [(0, 1), (1, 2)],
                                            fixations,
                                            ftimes, 2)

    def testPickle(self):
        p = pickle.dumps(self.pop)
        up = pickle.loads(p)
        self.assertEqual(self.pop, up)


if __name__ == "__main__":
    unittest.main()
