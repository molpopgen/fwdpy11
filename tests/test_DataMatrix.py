import unittest
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
        self.hm_neutral = np.array(self.hm.neutral)
        self.hm_selected = np.array(self.hm.selected)
        self.gm_neutral = np.array(self.gm.neutral)
        self.gm_selected = np.array(self.gm.selected)

    def testNcol(self):
        self.assertEqual(self.hm.ncol, 2*len(self.indlist))
        self.assertEqual(self.gm.ncol, len(self.indlist))

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
        # Get the col sums
        colSums = self.hm_neutral.sum(axis=0)
        colSumsSel = self.hm_selected.sum(axis=0)
        self.assertEqual(len(colSums), self.hm_neutral.shape[1])
        self.assertEqual(len(colSumsSel), self.hm_selected.shape[1])
        j = 0
        for i in range(100, 150):
            # Num neutral variants in diploid i, gamete 0
            nmuts = len(self.pop.gametes[self.pop.diploids[i].first].mutations)
            self.assertEqual(nmuts, colSums[j])
            # Num neutral variants in diploid i, gamete 1
            nmuts = len(
                self.pop.gametes[self.pop.diploids[i].second].mutations)
            self.assertEqual(nmuts, colSums[j + 1])

            # Now, test numbers of selected
            nmuts = len(
                self.pop.gametes[self.pop.diploids[i].first].smutations)
            self.assertEqual(nmuts, colSumsSel[j])
            nmuts = len(
                self.pop.gametes[self.pop.diploids[i].second].smutations)
            self.assertEqual(nmuts, colSumsSel[j + 1])

            j += 2

    def testGenoMat(self):
        self.assertTrue(self.gm_neutral.dtype == np.int8)
        self.assertEqual(self.gm_neutral.ndim, 2)
        # Get the col sums
        colSums = self.gm_neutral.sum(axis=0)
        colSumsSel = self.gm_selected.sum(axis=0)
        self.assertEqual(len(colSums), self.gm_neutral.shape[1])
        self.assertEqual(len(colSumsSel), self.gm_selected.shape[1])
        j = 0
        for i in range(100, 150):
            nmuts = len(self.pop.gametes[self.pop.diploids[i].first].mutations)
            nmuts += len(self.pop.gametes[self.pop.diploids[i].second].mutations)
            self.assertEqual(nmuts, colSums[j])

            # Now, test numbers of selected
            nmuts = len(
                self.pop.gametes[self.pop.diploids[i].first].smutations)
            nmuts += len(self.pop.gametes[self.pop.diploids[i].second].smutations)
            self.assertEqual(nmuts, colSumsSel[j])

            j += 1

    def testHapMatToSample(self):
        neut_sample, sel_sample = fwdpy11.sampling.matrix_to_sample(self.hm)
        rowSums = self.hm_neutral.sum(axis=1)
        rowSumsSel = self.hm_selected.sum(axis=1)
        self.assertEqual(len(neut_sample), self.hm_neutral.shape[0])
        self.assertEqual(len(rowSums), len(neut_sample))
        self.assertEqual(len(sel_sample), self.hm_selected.shape[0])
        self.assertEqual(len(rowSumsSel), len(sel_sample))
        i = 0
        for ni in neut_sample:
            num_ones = ni[1].count('1')
            self.assertEqual(num_ones, rowSums[i])
            i += 1


class test_DataMatrixFromMlocusPop(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.pop = quick_mlocus_qtrait()
        self.indlist = [i for i in range(100, 150)]
        self.keys = fwdpy11.sampling.mutation_keys(self.pop, self.indlist)
        self.nkeys = sorted(
            self.keys[0], key=lambda x: self.pop.mutations[x[0]].pos)
        self.skeys = sorted(
            self.keys[1], key=lambda x: self.pop.mutations[x[0]].pos)
        self.hm = fwdpy11.sampling.haplotype_matrix(
            self.pop, self.indlist, self.nkeys, self.skeys)
        self.gm = fwdpy11.sampling.genotype_matrix(
            self.pop, self.indlist, self.nkeys, self.skeys)
        self.hm_neutral = np.array(self.hm.neutral)
        self.hm_selected = np.array(self.hm.selected)
        self.gm_neutral = np.array(self.gm.neutral)
        self.gm_selected = np.array(self.gm.selected)

    def testConvertHapMatrixToSample(self):
        nsample, ssample = fwdpy11.sampling.matrix_to_sample(self.hm)
        nsample_split = fwdpy11.sampling.separate_samples_by_loci(
            self.pop.locus_boundaries, nsample)
        self.assertEqual(len(nsample_split), len(self.pop.locus_boundaries))
        self.assertEqual(len(nsample_split), self.pop.nloci)
        self.assertEqual(len(nsample_split), len(self.pop.locus_boundaries))
        key = 0
        for locus_data in nsample_split:
            self.assertTrue(sorted(locus_data, key=lambda x: x[0]))
            for site in locus_data:
                self.assertTrue(site[0] >= self.pop.locus_boundaries[key][0])
                self.assertTrue(site[0] < self.pop.locus_boundaries[key][1])
                self.assertTrue(len(site[1]) == self.hm_neutral.shape[1])
            key += 1

    def testConvertGenoMatrixToSample(self):
        nsample, ssample = fwdpy11.sampling.matrix_to_sample(self.gm)
        nsample_split = fwdpy11.sampling.separate_samples_by_loci(
            self.pop.locus_boundaries, nsample)
        self.assertEqual(len(nsample_split), len(self.pop.locus_boundaries))
        key = 0
        for locus in nsample_split:
            for site in locus:
                self.assertTrue(site[0] >= self.pop.locus_boundaries[key][0])
                self.assertTrue(site[0] < self.pop.locus_boundaries[key][1])
                self.assertTrue(len(site[1]) == self.gm_neutral.shape[1])
            key += 1

    def testSampleByLocus(self):
        lm = self.pop.sample_by_locus(self.indlist)
        self.assertEqual(len(lm), self.pop.nloci)
        dm = self.pop.sample(self.indlist)
        nsample, ssample = fwdpy11.sampling.matrix_to_sample(dm)
        nsample_split = fwdpy11.sampling.separate_samples_by_loci(
            self.pop.locus_boundaries, nsample)
        for i, j in zip(lm, nsample_split):
            # Test that positons are the same
            pi = [k for k in i.neutral.positions]
            pj = [k[0] for k in j]
            self.assertEqual(pi, pj)
            im = np.array(i.neutral)

            # convert sample data to matrix-like data
            temp = []
            for k in j:
                for c in k[1]:
                    if c == '0':
                        temp.append(0)
                    else:
                        temp.append(1)

            tempa = np.array(temp, dtype=im.dtype).reshape(im.shape)
            self.assertTrue(np.array_equal(im, tempa))


if __name__ == "__main__":
    unittest.main()
