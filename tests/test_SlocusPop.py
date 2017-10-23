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


class testSlocusPop(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.pop = fp11.SlocusPop(1000)

    def test_N(self):
        self.assertEqual(self.pop.N, 1000)

    def test_generation(self):
        self.assertEqual(self.pop.generation, 0)

    def test_fitnesses(self):
        # All fitnesses should be 1
        self.assertEqual(
            sum([i.w for i in self.pop.diploids]),
            float(self.pop.N))

    def test_labels(self):
        # All diploids should be labeled 0 to pop.N-1
        self.assertEqual(
            [i.label for i in self.pop.diploids],
            [i for i in range(self.pop.N)])

    def test_genetic_values(self):
        self.assertEqual(sum([i.g for i in self.pop.diploids]), 0.0)

    def test_e_values(self):
        self.assertEqual(sum([i.e for i in self.pop.diploids]), 0.0)


class testSlocusPopExceptions(unittest.TestCase):
    def testNzero(self):
        with self.assertRaises(ValueError):
            fp11.SlocusPop(0)


class testSampling(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        from quick_pops import quick_nonneutral_slocus
        self.pop = quick_nonneutral_slocus()
        self.rng = fp11.GSLrng(42)

    def testRandomSample(self):
        self.pop.sample(self.rng, 10)
        self.pop.sample(self.rng, 10, separate=False)
        self.pop.sample(self.rng, 10, remove_fixed=False)
        self.pop.sample(self.rng, 10, separate=True, remove_fixed=False)
        self.pop.sample(self.rng, 10, separate=False, remove_fixed=False)
        self.pop.sample(self.rng, 10, separate=True, remove_fixed=True)

    def testDefinedSample(self):
        self.pop.sample_ind(range(10))

        with self.assertRaises(IndexError):
            """
            fwdpp catches case where i >= N
            """
            self.pop.sample_ind(range(self.pop.N, self.pop.N + 10))

        with self.assertRaises(TypeError):
            """
            pybind11 disallows conversion of negative
            numbers to a list of unsigned types.
            """
            self.pop.sample_ind(range(-10, 10))


if __name__ == "__main__":
    unittest.main()
