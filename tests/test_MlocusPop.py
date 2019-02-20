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
import fwdpy11 as fp11


class testMlocusPop(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.pop = fp11.MlocusPop(1000, [(i, i + 1) for i in range(5)])

    def test_N(self):
        self.assertEqual(self.pop.N, 1000)

    def test_nloci(self):
        self.assertEqual(self.pop.nloci, 5)

    def test_fitnesses(self):
        # All fitnesses should be 1
        self.assertEqual(
            sum([i.w for i in self.pop.diploid_metadata]),
            float(self.pop.N))

    def test_labels(self):
        # All diploids should be labeled 0 to pop.N-1
        self.assertEqual(
            [i.label for i in self.pop.diploid_metadata],
            [i for i in range(self.pop.N)])

    def test_genetic_values(self):
        self.assertEqual(sum([i.g for i in self.pop.diploid_metadata]), 0.0)

    def test_e_values(self):
        self.assertEqual(sum([i.e for i in self.pop.diploid_metadata]), 0.0)

    def test_nodes(self):
        self.assertEqual(len(self.pop.diploids),
                         len(self.pop.diploid_metadata))
        import fwdpy11.ts
        for i in self.pop.diploid_metadata:
            self.assertEqual(i.nodes[0], fwdpy11.ts.NULL_NODE)
            self.assertEqual(i.nodes[1], fwdpy11.ts.NULL_NODE)


class testMlocusPopExceptions(unittest.TestCase):
    def testNzero(self):
        with self.assertRaises(ValueError):
            fp11.MlocusPop(0, [(i, i + 1) for i in range(5)])

    def testNoLoci(self):
        with self.assertRaises(ValueError):
            fp11.MlocusPop(1000, [])


class testSampling(unittest.TestCase):
    @classmethod
    def setUp(self):
        from quick_pops import quick_mlocus_qtrait
        self.pop = quick_mlocus_qtrait()
        self.rng = fp11.GSLrng(42)

    def testRandomSample(self):
        x = self.pop.sample(rng=self.rng, nsam=10)
        # self.assertEqual(len(x), self.pop.nloci)
        # x = self.pop.sample(rng=self.rng, nsam=10, separate=False)
        # self.assertEqual(len(x), self.pop.nloci)
        # x = self.pop.sample(rng=self.rng, nsam=10, remove_fixed=False)
        # self.assertEqual(len(x), self.pop.nloci)
        # x = self.pop.sample(rng=self.rng, nsam=10,
        #                     separate=True, remove_fixed=False)
        # self.assertEqual(len(x), self.pop.nloci)
        # for i in x:
        #     self.assertTrue(isinstance(i, tuple))
        #     self.assertEqual(len(i), 2)
        # x = self.pop.sample(rng=self.rng, nsam=10,
        #                     separate=False, remove_fixed=False)
        # self.assertEqual(len(x), self.pop.nloci)
        # x = self.pop.sample(rng=self.rng, nsam=10,
        #                     separate=True, remove_fixed=True)
        # self.assertEqual(len(x), self.pop.nloci)

    def testDefinedSample(self):
        x = self.pop.sample(individuals=range(10))
        #self.assertEqual(len(x), self.pop.nloci)
        with self.assertRaises(IndexError):
            """
            fwdpp catches case where i >= N
            """
            self.pop.sample(individuals=range(self.pop.N, self.pop.N + 10))

        with self.assertRaises(Exception):
            """
            pybind11 disallows conversion of negative
            numbers to a list of unsigned types.
            """
            self.pop.sample(individuals=range(-10, 10))

    def test_nodes(self):
        import fwdpy11.ts
        for i in self.pop.diploid_metadata:
            self.assertEqual(i.nodes[0], fwdpy11.ts.NULL_NODE)
            self.assertEqual(i.nodes[1], fwdpy11.ts.NULL_NODE)


class testBadLocusBoundaries(unittest.TestCase):
    def testUnsorted(self):
        with self.assertRaises(ValueError):
            pop = fp11.MlocusPop(100, [(1, 2), (0, 1)])

    def testBadInterval(self):
        with self.assertRaises(ValueError):
            pop = fp11.MlocusPop(100, [(0, 1), (1, 1)])
        with self.assertRaises(ValueError):
            pop = fp11.MlocusPop(100, [(0, 1), (2, 1)])

    def testOverlap(self):
        with self.assertRaises(ValueError):
            pop = fp11.MlocusPop(100, [(0, 1.1), (1, 2)])
        with self.assertRaises(ValueError):
            pop = fp11.MlocusPop(100, [(0, 3), (1, 2)])


class testSettingLocusBoundaries(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.pop = fp11.MlocusPop(100, [(0, 1), (1, 2)])

    def testReassign(self):
        self.pop.locus_boundaries = [(2, 3), (11, 45)]

    def testBadReassign(self):
        with self.assertRaises(ValueError):
            self.pop.locus_boundaries = []
        with self.assertRaises(ValueError):
            self.pop.locus_boundaries = [(1, 2.), (1, 3)]


if __name__ == "__main__":
    unittest.main()
