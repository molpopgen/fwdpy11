import unittest
import fwdpy11 as fp11
import fwdpy11.sampling
import numpy as np
import numpy as np
import pyximport
pyximport.install(setup_args={'include_dirs': np.get_include()})
from mean_sel_coeff import mean_s


class test_SlocusPop(unittest.TestCase):
    """
    Test coercion of population data
    to numpy arrays
    """
    @classmethod
    def setUpClass(self):
        from quick_pops import quick_nonneutral_slocus
        self.pop = quick_nonneutral_slocus()
        self.muts = np.array(self.pop.mutations.array())

    def testMutations(self):
        self.assertTrue(isinstance(self.muts, np.ndarray))
        self.assertEqual(self.muts.shape[0], len(self.pop.mutations))
        for i, j in zip(self.muts, self.pop.mutations):
            self.assertEqual(i['neutral'], j.neutral)
            self.assertEqual(i['h'], j.h)
            self.assertEqual(i['pos'], j.pos)
            self.assertEqual(i['g'], j.g)
            self.assertEqual(i['s'], j.s)

    def testCythonFunc(self):
        ms = mean_s(self.muts['s'])
        self.assertAlmostEqual(ms, self.muts['s'].sum())

    def testDiploidTraits(self):
        dips = np.array(self.pop.diploids.trait_array())
        for i, j in zip(dips, self.pop.diploids):
            self.assertEqual(i['g'], j.g)
            self.assertEqual(i['e'], j.e)
            self.assertEqual(i['w'], j.w)

    def testDiploidTraitsSlice(self):
        dips = np.array(self.pop.diploids.trait_array(slice(15, 103, 11)))
        for i, j in zip(dips, self.pop.diploids[15:103:11]):
            self.assertEqual(i['g'], j.g)
            self.assertEqual(i['e'], j.e)
            self.assertEqual(i['w'], j.w)

    def testDiploidKeys(self):
        keys = np.array(self.pop.diploids.key_array())
        for i, j in zip(keys, self.pop.diploids):
            self.assertEqual(i['first'], j.first)
            self.assertEqual(i['second'], j.second)


class test_MlocusPop(unittest.TestCase):
    """
    Test coercion of population data
    to numpy arrays
    """
    @classmethod
    def setUpClass(self):
        from quick_pops import quick_mlocus_qtrait
        self.pop = quick_mlocus_qtrait()
        self.muts = np.array(self.pop.mutations.array())

    def testMutations(self):
        self.assertTrue(isinstance(self.muts, np.ndarray))
        self.assertEqual(self.muts.shape[0], len(self.pop.mutations))
        for i, j in zip(self.muts, self.pop.mutations):
            self.assertEqual(i['neutral'], j.neutral)
            self.assertEqual(i['h'], j.h)
            self.assertEqual(i['pos'], j.pos)
            self.assertEqual(i['g'], j.g)
            self.assertEqual(i['s'], j.s)

    def testDiploidTraits(self):
        dips = np.array(self.pop.diploids.trait_array())
        for i, j in zip(dips, self.pop.diploids):
            self.assertEqual(i['g'], j[0].g)
            self.assertEqual(i['e'], j[0].e)
            self.assertEqual(i['w'], j[0].w)

    def testDiploidTraitsSlice(self):
        dips = np.array(self.pop.diploids.trait_array(slice(15, 103, 11)))
        for i, j in zip(dips, self.pop.diploids[15:103:11]):
            self.assertEqual(i['g'], j[0].g)
            self.assertEqual(i['e'], j[0].e)
            self.assertEqual(i['w'], j[0].w)

    def testDiploidKeys(self):
        keys = np.array(self.pop.diploids.key_array())
        keys_per_dip_offsets = [i for i in range(len(keys))[::self.pop.nloci]]
        self.assertEqual(len(keys_per_dip_offsets), len(self.pop.diploids))
        for i, j in zip(keys_per_dip_offsets, self.pop.diploids):
            for k, l in zip(keys[i:i + self.pop.nloci], j):
                self.assertEqual(k['first'], l.first)
                self.assertEqual(k['second'], l.second)

if __name__ == "__main__":
    unittest.main()

