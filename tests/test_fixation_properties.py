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
# Test behavior of fixations in simulations.
# As of 0.1.3a1, fixations are to be sorted by position
# for all types of simulations, allowing for easier lookup.
#
# In 0.1.5, this unit test was refactored to skip simulation.
# Instead, we use pybind11 + cppimport to expose the underlying
# C++ API to Python, create pops manually, and test the API.

import unittest
import os
import fwdpy11
import cppimport
cppimport.force_rebuild()
fp = cppimport.imp("fixation_properties")


class testFixationsAreSortedSlocusPop(unittest.TestCase):
    @classmethod
    def setUp(self):
        mutations = fwdpy11.VecMutation()
        fixations = fwdpy11.VecMutation()
        gametes = fwdpy11.VecGamete()
        diploids = fwdpy11.VecDiploid()
        mutations.append(fwdpy11.Mutation(0.1, -0.01, 1.0, 0, 0))
        mutations.append(fwdpy11.Mutation(0.2, 0.0, 1.0, 0, 0))
        gametes.append(fwdpy11.Gamete(
            (4, fwdpy11.VecUint32([1]), fwdpy11.VecUint32([0]))))
        diploids.append(fwdpy11.SingleLocusDiploid(0, 0))
        diploids.append(fwdpy11.SingleLocusDiploid(0, 0))
        self.pop = fwdpy11.SlocusPop(diploids, gametes, mutations)

    def testSetup(self):
        self.assertEqual(len(self.pop.mutations), 2)
        self.assertEqual(len(self.pop.gametes), 1)
        self.assertEqual(len(self.pop.diploids), 2)
        for g in self.pop.gametes:
            self.assertEqual(len(g.mutations), 1)
            self.assertEqual(len(g.smutations), 1)
        for mc in self.pop.mcounts:
            self.assertEqual(mc, 4)

    def testWithPruneSelected(self):
        fp.update_mutations(self.pop, 1, 2 * len(self.pop.diploids), True)
        self.assertEqual(len(self.pop.fixations), 2)
        self.assertEqual(len(self.pop.fixation_times), 2)
        self.assertEqual(len([i for i in self.pop.mcounts if i == 0]), 2)
        self.assertTrue(sorted(self.pop.fixations, key=lambda x: x.pos))
        for mc in self.pop.mcounts:
            self.assertEqual(mc, 0)

    def testWithoutPruneSelected(self):
        fp.update_mutations(self.pop, 1, 2 * len(self.pop.diploids), False)
        self.assertEqual(len(self.pop.fixations), 2)
        self.assertEqual(len(self.pop.fixation_times), 2)
        self.assertEqual(len([i for i in self.pop.mcounts if i == 0]), 1)
        self.assertTrue(sorted(self.pop.fixations, key=lambda x: x.pos))
        for m, mc in zip(self.pop.mutations, self.pop.mcounts):
            if m.neutral is True:
                self.assertEqual(mc, 0)
            else:
                self.assertEqual(mc, 4)


class testFixationsAreSortedMlocusPop(unittest.TestCase):
    @classmethod
    def setUp(self):
        mutations = fwdpy11.VecMutation()
        fixations = fwdpy11.VecMutation()
        gametes = fwdpy11.VecGamete()
        diploids = fwdpy11.VecVecDiploid()
        mutations.append(fwdpy11.Mutation(0.25, -0.01, 1.0, 0, 0))
        mutations.append(fwdpy11.Mutation(0.75, 0.0, 1.0, 0, 0))
        mutations.append(fwdpy11.Mutation(1.25, 0.0, 1.0, 0, 0))
        mutations.append(fwdpy11.Mutation(1.75, 0.25, 1.0, 0, 0))
        gametes.append(fwdpy11.Gamete(
            (4, fwdpy11.VecUint32([1]), fwdpy11.VecUint32([0]))))
        gametes.append(fwdpy11.Gamete(
            (4, fwdpy11.VecUint32([2]), fwdpy11.VecUint32([3]))))
        dip = fwdpy11.VecDiploid(
            [fwdpy11.SingleLocusDiploid(0, 0), fwdpy11.SingleLocusDiploid(1, 1)])
        diploids.append(dip)
        diploids.append(dip)
        self.pop = fwdpy11.MlocusPop(
            diploids, gametes, mutations, [(0, 1), (1, 2)])

    def testSetup(self):
        self.assertEqual(len(self.pop.mutations), 4)
        self.assertEqual(len(self.pop.gametes), 2)
        self.assertEqual(len(self.pop.diploids), 2)
        for d in self.pop.diploids:
            self.assertEqual(len(d), 2)
            for l in d:
                self.assertEqual(len(self.pop.gametes[l.first].mutations), 1)
                self.assertEqual(len(self.pop.gametes[l.first].smutations), 1)
                self.assertEqual(len(self.pop.gametes[l.second].mutations), 1)
                self.assertEqual(len(self.pop.gametes[l.second].smutations), 1)

        self.assertEqual(len(self.pop.locus_boundaries), 2)
        self.assertEqual(self.pop.nloci, 2)
        for g in self.pop.gametes:
            self.assertEqual(len(g.mutations), 1)
            self.assertEqual(len(g.smutations), 1)
        for mc in self.pop.mcounts:
            self.assertEqual(mc, 4)

    def testWithPruneSelected(self):
        fp.update_mutations(self.pop, 1, 2 * len(self.pop.diploids), True)
        self.assertEqual(len(self.pop.fixations), 4)
        self.assertEqual(len(self.pop.fixation_times), 4)
        self.assertEqual(len([i for i in self.pop.mcounts if i == 0]), 4)
        self.assertTrue(sorted(self.pop.fixations, key=lambda x: x.pos))
        for mc in self.pop.mcounts:
            self.assertEqual(mc, 0)

    def testWithoutPruneSelected(self):
        fp.update_mutations(self.pop, 1, 2 * len(self.pop.diploids), False)
        self.assertEqual(len(self.pop.fixations), 4)
        self.assertEqual(len(self.pop.fixation_times), 4)
        self.assertEqual(len([i for i in self.pop.mcounts if i == 0]), 2)
        self.assertTrue(sorted(self.pop.fixations, key=lambda x: x.pos))
        for m, mc in zip(self.pop.mutations, self.pop.mcounts):
            if m.neutral is True:
                self.assertEqual(mc, 0)
            else:
                self.assertEqual(mc, 4)


if __name__ == "__main__":
    unittest.main()
