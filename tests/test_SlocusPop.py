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
        x = self.pop.sample(rng=self.rng, nsam=10)
        self.assertTrue(type(x) is tuple)
        x = self.pop.sample(rng=self.rng, nsam=10, separate=False)
        self.assertTrue(type(x) is list)
        x = self.pop.sample(rng=self.rng, nsam=10, remove_fixed=False)
        self.assertTrue(type(x) is tuple)
        x = self.pop.sample(rng=self.rng, nsam=10,
                            separate=True, remove_fixed=False)
        self.assertTrue(type(x) is tuple)

    def testDefinedSample(self):
        self.pop.sample(individuals=range(10))

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


class testPythonObjects(unittest.TestCase):
    @classmethod
    def setUp(self):
        from quick_pops import quick_slocus_qtrait_pop_params
        self.pop, self.pdict = quick_slocus_qtrait_pop_params()
        self.rng = fp11.GSLrng(101)

    def testInitialState(self):
        self.assertTrue(self.pop.popdata is None)
        self.assertTrue(self.pop.popdata_user is None)

    def testParentalData(self):
        from fwdpy11.model_params import SlocusParamsQ
        from fwdpy11.wright_fisher_qtrait import evolve
        params = SlocusParamsQ(**self.pdict)
        evolve(self.rng, self.pop, params)
        parents = [i.parental_data for i in self.pop.diploids]
        for i in parents:
            self.assertTrue(i is not None)
            self.assertTrue(len(i) == 2)
            self.assertTrue(i[0] < self.pop.N)
            self.assertTrue(i[1] < self.pop.N)

    def testPopdataReadOnly(self):
        from fwdpy11.model_params import SlocusParamsQ
        from fwdpy11.wright_fisher_qtrait import evolve
        params = SlocusParamsQ(**self.pdict)
        with self.assertRaises(Exception):
            class Recorder(object):
                def __call__(self, pop):
                    pop.popdata = []
            r = Recorder()
            evolve(self.rng, self.pop, params, r)

    def testPopdataUserAndPickling(self):
        from fwdpy11.model_params import SlocusParamsQ
        from fwdpy11.wright_fisher_qtrait import evolve
        params = SlocusParamsQ(**self.pdict)

        class Recorder(object):
            def __call__(self, pop):
                pop.popdata_user = pop.generation
        r = Recorder()
        evolve(self.rng, self.pop, params, r)
        self.assertEqual(self.pop.generation, self.pop.popdata_user)

        import pickle
        d = pickle.dumps(self.pop)
        up = pickle.loads(d)
        self.assertEqual(up.popdata_user, self.pop.generation)
        self.assertEqual(self.pop, up)


if __name__ == "__main__":
    unittest.main()
