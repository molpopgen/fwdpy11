import unittest
import fwdpy11
import numpy as np
import pyximport
pyximport.install(setup_args={'include_dirs': np.get_include()})
from mean_sel_coeff import mean_s


class test_DiploidPopulation(unittest.TestCase):
    """
    Test coercion of population data
    to numpy arrays
    """
    @classmethod
    def setUpClass(self):
        from quick_pops import quick_nonneutral_slocus
        self.pop = quick_nonneutral_slocus()
        self.muts = self.pop.mutations_ndarray

    def testMutations(self):
        self.assertTrue(isinstance(self.muts, np.ndarray))
        self.assertEqual(self.muts.shape[0], len(self.pop.mutations))
        for i, j in zip(self.muts, self.pop.mutations):
            self.assertEqual(i['neutral'], j.neutral)
            self.assertEqual(i['h'], j.h)
            self.assertEqual(i['pos'], j.pos)
            self.assertEqual(i['g'], j.g)
            self.assertEqual(i['s'], j.s)

            m = fwdpy11.Mutation(tuple(i)[:-1])
            self.assertEqual(m.pos, j.pos)
            self.assertEqual(m.s, j.s)
            self.assertEqual(m.h, j.h)
            self.assertEqual(m.g, j.g)
            self.assertEqual(m.label, j.label)
            self.assertEqual(m.neutral, j.neutral)

    def testCythonFunc(self):
        ms = mean_s(self.muts['s'])
        self.assertAlmostEqual(ms, self.muts['s'].mean())

    def testDiploidTraits(self):
        dips = np.array(self.pop.diploid_metadata)
        for i, j in zip(dips, self.pop.diploid_metadata):
            self.assertEqual(i['g'], j.g)
            self.assertEqual(i['e'], j.e)
            self.assertEqual(i['w'], j.w)

    def testDiploidTraitsSlice(self):
        dips = np.array(self.pop.diploid_metadata[15: 103: 11])
        for i, j in zip(dips, self.pop.diploid_metadata[15:103:11]):
            self.assertEqual(i['g'], j.g)
            self.assertEqual(i['e'], j.e)
            self.assertEqual(i['w'], j.w)

    def testDiploidKeys(self):
        keys = np.array(self.pop.diploids)
        for i, j in zip(keys, self.pop.diploids):
            self.assertEqual(i['first'], j.first)
            self.assertEqual(i['second'], j.second)


if __name__ == "__main__":
    unittest.main()
