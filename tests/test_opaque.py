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

# This file tests that types in population objects
# are indeed "opaque" and not being copied to/from
# Python lists

import unittest
import fwdpy11


class testSlocusPop(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.pop = fwdpy11.SlocusPop(1000)

    def test_mutations(self):
        self.assertTrue(
            type(self.pop.mutations)
            is fwdpy11.fwdpy11_types.MutationContainer)
        self.assertFalse(type(self.pop.mutations) is list)
        self.assertTrue(type(list(self.pop.mutations)) is list)

    def test_mcount(self):
        self.assertTrue(
            type(self.pop.mcounts)
            is fwdpy11.fwdpy11_types.VectorUint32)

    def test_diploids(self):
        self.assertTrue(
            type(self.pop.diploids)
            is fwdpy11.fwdpy11_types.DiploidContainer)

    def test_a_diploid(self):
        self.assertTrue(
            type(self.pop.diploids[0])
            is fwdpy11.fwdpy11_types.SingleLocusDiploid)

    def test_gametes(self):
        self.assertTrue(
            type(self.pop.gametes)
            is fwdpy11.fwdpy11_types.GameteContainer)


class SlocusTypeSampler(object):
    def __call__(self, pop):
        assert(type(pop.mutations) is fwdpy11.fwdpy11_types.MutationContainer)
        assert(type(pop.mutations) is not list)
        assert(type(list(pop.mutations)) is list)
        assert(type(pop.mcounts) is fwdpy11.fwdpy11_types.VectorUint32)
        assert(type(pop.diploids) is fwdpy11.fwdpy11_types.DiploidContainer)
        assert(type(pop.diploids[0])
               is fwdpy11.fwdpy11_types.SingleLocusDiploid)
        assert(type(pop.gametes) is fwdpy11.fwdpy11_types.GameteContainer)


class testSlocusPopSampler(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        from fwdpy11.ezparams import mslike
        from fwdpy11.model_params import SlocusParams
        self.sampler = SlocusTypeSampler()
        self.pop = fwdpy11.SlocusPop(1000)
        self.params_dict = mslike(self.pop, simlen=10)
        self.params = SlocusParams(**self.params_dict)
        self.rng = fwdpy11.GSLrng(42)

    def testSampler(self):
        from fwdpy11.wright_fisher import evolve
        try:
            evolve(self.rng, self.pop, self.params, self.sampler)
        except:
            self.fail("unexpcted AssertionError")


class testMlocusPop(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.pop = fwdpy11.MlocusPop(1000, 5)

    def test_mutations(self):
        self.assertTrue(
            type(self.pop.mutations)
            is fwdpy11.fwdpy11_types.MutationContainer)
        self.assertFalse(type(self.pop.mutations) is list)
        self.assertTrue(type(list(self.pop.mutations)) is list)

    def test_mcount(self):
        self.assertTrue(
            type(self.pop.mcounts)
            is fwdpy11.fwdpy11_types.VectorUint32)

    def test_diploids(self):
        self.assertFalse(type(self.pop.diploids) is list)
        self.assertTrue(
            type(self.pop.diploids)
            is fwdpy11.fwdpy11_types.VecDiploidContainer)

    def test_a_diploid(self):
        self.assertTrue(
            type(self.pop.diploids[0])
            is fwdpy11.fwdpy11_types.DiploidContainer)

    def test_gametes(self):
        self.assertTrue(
            type(self.pop.gametes)
            is fwdpy11.fwdpy11_types.GameteContainer)

if __name__ == "__main__":
    unittest.main()

