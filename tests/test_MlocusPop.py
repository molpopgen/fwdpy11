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
        self.pop = fp11.MlocusPop(1000, 5)

    def test_N(self):
        self.assertEqual(self.pop.N, 1000)

    def test_nloci(self):
        self.assertEqual(self.pop.nloci, 5)

    def test_fitnesses(self):
        # All fitnesses should be 1
        self.assertEqual(
            sum([i[0].w for i in self.pop.diploids]),
            float(self.pop.N))

    def test_labels(self):
        # All diploids should be labeled 0 to pop.N-1
        self.assertEqual(
            [i[0].label for i in self.pop.diploids],
            [i for i in range(self.pop.N)])

    def test_genetic_values(self):
        self.assertEqual(sum([i[0].g for i in self.pop.diploids]), 0.0)

    def test_e_values(self):
        self.assertEqual(sum([i[0].e for i in self.pop.diploids]), 0.0)


class testMlocusPopExceptions(unittest.TestCase):
    def testNzero(self):
        with self.assertRaises(ValueError):
            p = fp11.MlocusPop(0, 5)

    def testNoLoci(self):
        with self.assertRaises(ValueError):
            p = fp11.MlocusPop(1000, 0)

if __name__ == "__main__":
    unittest.main()
