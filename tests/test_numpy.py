import unittest

import numpy as np

import fwdpy11


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
            self.assertEqual(i["neutral"], j.neutral)
            self.assertEqual(i["h"], j.h)
            self.assertEqual(i["pos"], j.pos)
            self.assertEqual(i["g"], j.g)
            self.assertEqual(i["s"], j.s)

    def testDiploidTraits(self):
        dips = np.array(self.pop.diploid_metadata)
        for i, j in zip(dips, self.pop.diploid_metadata):
            self.assertEqual(i["g"], j.g)
            self.assertEqual(i["e"], j.e)
            self.assertEqual(i["w"], j.w)

    def testDiploidTraitsSlice(self):
        dips = np.array(self.pop.diploid_metadata[15:103:11])
        for i, j in zip(dips, self.pop.diploid_metadata[15:103:11]):
            self.assertEqual(i["g"], j.g)
            self.assertEqual(i["e"], j.e)
            self.assertEqual(i["w"], j.w)

    def testDiploidKeys(self):
        keys = np.array(self.pop.diploids)
        for i, j in zip(keys, self.pop.diploids):
            self.assertEqual(i["first"], j.first)
            self.assertEqual(i["second"], j.second)


if __name__ == "__main__":
    unittest.main()
