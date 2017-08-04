import unittest
import fwdpy11 as fp11
import fwdpy11.sampling
import numpy as np
from quick_pops import quick_nonneutral_slocus, quick_mlocus_qtrait


class test_DataMatrixFromSlocusPop(unittest.TestCase):
    """
    Much of this is already tested in fwdpp's unit tests,
    but the purpose here is to make sure that the wrappers
    are working appropriately.
    """
    @classmethod
    def setUpClass(self):
        self.pop = quick_nonneutral_slocus()
        self.indlist = [i for i in range(100, 150)]
        self.keys = fwdpy11.sampling.mutation_keys(self.pop, self.indlist)
        self.hm = fwdpy11.sampling.haplotype_matrix(
            self.pop, self.indlist, self.keys[0], self.keys[1])
        self.gm = fwdpy11.sampling.genotype_matrix(
            self.pop, self.indlist, self.keys[0], self.keys[1])
        self.hm_neutral = np.ndarray(
            self.hm.ndim_neutral(),
            buffer=self.hm.neutral, dtype=np.int8)
        self.hm_selected = np.ndarray(
            self.hm.ndim_selected(),
            buffer=self.hm.selected, dtype=np.int8)
        self.gm_neutral = np.ndarray(
            self.gm.ndim_neutral(),
            buffer=self.gm.neutral, dtype=np.int8)
        self.gm_selected = np.ndarray(
            self.gm.ndim_selected(),
            buffer=self.gm.selected, dtype=np.int8)

    def testKeyNeutralityAndCount(self):
        for i in self.keys[0]:
            self.assertTrue(self.pop.mutations[i[0]].neutral)
            self.assertTrue(self.pop.mcounts[i[0]] > 0)
        for i in self.keys[1]:
            self.assertFalse(self.pop.mutations[i[0]].neutral)
            self.assertTrue(self.pop.mcounts[i[0]] > 0)

    def testHapMat(self):
        self.assertTrue(self.hm_neutral.dtype == np.int8)
        self.assertEqual(self.hm_neutral.ndim, 2)
        self.assertTrue(self.hm_selected.dtype == np.int8)
        self.assertEqual(self.hm_selected.ndim, 2)
        # Get the row sums
        rowSums = self.hm_neutral.sum(axis=1)
        rowSumsSel = self.hm_selected.sum(axis=1)
        self.assertEqual(len(rowSums), self.hm_neutral.shape[0])
        self.assertEqual(len(rowSumsSel), self.hm_selected.shape[0])
        j = 0
        for i in range(100, 150):
            # Num neutral variants in diploid i, gamete 0
            nmuts = len(self.pop.gametes[self.pop.diploids[i].first].mutations)
            self.assertEqual(nmuts, rowSums[j])
            # Num neutral variants in diploid i, gamete 1
            nmuts = len(
                self.pop.gametes[self.pop.diploids[i].second].mutations)
            self.assertEqual(nmuts, rowSums[j + 1])

            # Now, test numbers of selected
            nmuts = len(
                self.pop.gametes[self.pop.diploids[i].first].smutations)
            self.assertEqual(nmuts, rowSumsSel[j])
            nmuts = len(
                self.pop.gametes[self.pop.diploids[i].second].smutations)
            self.assertEqual(nmuts, rowSumsSel[j + 1])

            j += 2

    def testGenoMat(self):
        self.assertTrue(self.gm_neutral.dtype == np.int8)
        self.assertEqual(self.gm_neutral.ndim, 2)
        # Get the row sums
        rowSums = self.gm_neutral.sum(axis=1)
        rowSumsSel = self.gm_selected.sum(axis=1)
        self.assertEqual(len(rowSums), self.gm_neutral.shape[0])
        self.assertEqual(len(rowSumsSel), self.gm_selected.shape[0])
        j = 0
        for i in range(100, 150):
            nmuts = len(self.pop.gametes[self.pop.diploids[i].first].mutations)
            nmuts += len(self.pop.gametes[self.pop.diploids[i].second].mutations)
            self.assertEqual(nmuts, rowSums[j])

            # Now, test numbers of selected
            nmuts = len(
                self.pop.gametes[self.pop.diploids[i].first].smutations)
            nmuts += len(self.pop.gametes[self.pop.diploids[i].second].smutations)
            self.assertEqual(nmuts, rowSumsSel[j])

            j += 1

    def testHapMatToSample(self):
        neut_sample = fwdpy11.sampling.matrix_to_sample(self.hm, True)
        sel_sample = fwdpy11.sampling.matrix_to_sample(self.hm, False)
        colSums = self.hm_neutral.sum(axis=0)
        colSumsSel = self.hm_selected.sum(axis=0)
        self.assertEqual(len(neut_sample), self.hm_neutral.shape[1])
        self.assertEqual(len(colSums), len(neut_sample))
        self.assertEqual(len(sel_sample), self.hm_selected.shape[1])
        self.assertEqual(len(colSumsSel), len(sel_sample))
        i = 0
        for ni in neut_sample:
            num_ones = ni[1].count('1')
            self.assertEqual(num_ones, colSums[i])
            i += 1


class test_DataMatrixFromMlocusPop(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.pop = quick_mlocus_qtrait()
        self.indlist = [i for i in range(100, 150)]
        self.keys = fwdpy11.sampling.mutation_keys(self.pop, self.indlist)
        self.hm = fwdpy11.sampling.haplotype_matrix(
            self.pop, self.indlist, self.keys[0], self.keys[1])
        self.gm = fwdpy11.sampling.genotype_matrix(
            self.pop, self.indlist, self.keys[0], self.keys[1])
        self.hm_neutral = np.ndarray(
            self.hm.ndim_neutral(),
            buffer=self.hm.neutral, dtype=np.int8)
        self.hm_selected = np.ndarray(
            self.hm.ndim_selected(),
            buffer=self.hm.selected, dtype=np.int8)
        self.gm_neutral = np.ndarray(
            self.gm.ndim_neutral(),
            buffer=self.gm.neutral, dtype=np.int8)
        self.gm_selected = np.ndarray(
            self.gm.ndim_selected(),
            buffer=self.gm.selected, dtype=np.int8)

    def testConvertHapMatrixToSample(self):
        nsample = fwdpy11.sampling.matrix_to_sample(self.hm, True)
        nsample_split = fwdpy11.sampling.separate_samples_by_loci(
            self.pop.locus_boundaries, nsample)
        self.assertEqual(len(nsample_split), len(self.pop.locus_boundaries))
        for key, value in nsample_split.items():
            self.assertTrue(key < self.pop.nloci)
            for site in value:
                self.assertTrue(site[0] >= self.pop.locus_boundaries[key][0])
                self.assertTrue(site[0] < self.pop.locus_boundaries[key][1])
                self.assertTrue(len(site[1]) == self.hm_neutral.shape[0])

    def testConvertGenoMatrixToSample(self):
        nsample = fwdpy11.sampling.matrix_to_sample(self.gm, True)
        nsample_split = fwdpy11.sampling.separate_samples_by_loci(
            self.pop.locus_boundaries, nsample)
        self.assertEqual(len(nsample_split), len(self.pop.locus_boundaries))
        for key, value in nsample_split.items():
            self.assertTrue(key < self.pop.nloci)
            for site in value:
                self.assertTrue(site[0] >= self.pop.locus_boundaries[key][0])
                self.assertTrue(site[0] < self.pop.locus_boundaries[key][1])
                self.assertTrue(len(site[1]) == self.gm_neutral.shape[0])

if __name__ == "__main__":
    unittest.main()
