import unittest

import numpy as np

import fwdpy11
from quick_pops import quick_nonneutral_slocus


class TestDataMatrixFromDiploidPopulation(unittest.TestCase):
    """
    Much of this is already tested in fwdpp's unit tests,
    but the purpose here is to make sure that the wrappers
    are working appropriately.
    """

    @classmethod
    def setUpClass(self):
        self.pop = quick_nonneutral_slocus()
        self.indlist = [i for i in range(100, 150)]
        self.keys = fwdpy11.mutation_keys(self.pop, self.indlist)
        self.hm = fwdpy11.haplotype_matrix(
            self.pop, self.indlist, self.keys[0], self.keys[1]
        )
        self.gm = fwdpy11.genotype_matrix(
            self.pop, self.indlist, self.keys[0], self.keys[1]
        )
        self.hm_neutral = np.array(self.hm.neutral)
        self.hm_selected = np.array(self.hm.selected)
        self.gm_neutral = np.array(self.gm.neutral)
        self.gm_selected = np.array(self.gm.selected)

    def testNcol(self):
        self.assertEqual(self.hm.ncol, 2 * len(self.indlist))
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
            nmuts = len(self.pop.haploid_genomes[self.pop.diploids[i].first].mutations)
            self.assertEqual(nmuts, colSums[j])
            # Num neutral variants in diploid i, gamete 1
            nmuts = len(self.pop.haploid_genomes[self.pop.diploids[i].second].mutations)
            self.assertEqual(nmuts, colSums[j + 1])

            # Now, test numbers of selected
            nmuts = len(self.pop.haploid_genomes[self.pop.diploids[i].first].smutations)
            self.assertEqual(nmuts, colSumsSel[j])
            nmuts = len(
                self.pop.haploid_genomes[self.pop.diploids[i].second].smutations
            )
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
            nmuts = len(self.pop.haploid_genomes[self.pop.diploids[i].first].mutations)
            nmuts += len(
                self.pop.haploid_genomes[self.pop.diploids[i].second].mutations
            )
            self.assertEqual(nmuts, colSums[j])

            # Now, test numbers of selected
            nmuts = len(self.pop.haploid_genomes[self.pop.diploids[i].first].smutations)
            nmuts += len(
                self.pop.haploid_genomes[self.pop.diploids[i].second].smutations
            )
            self.assertEqual(nmuts, colSumsSel[j])

            j += 1

    def testHapMatToSample(self):
        neut_sample, sel_sample = fwdpy11.matrix_to_sample(self.hm)
        rowSums = self.hm_neutral.sum(axis=1)
        rowSumsSel = self.hm_selected.sum(axis=1)
        self.assertEqual(len(neut_sample), self.hm_neutral.shape[0])
        self.assertEqual(len(rowSums), len(neut_sample))
        self.assertEqual(len(sel_sample), self.hm_selected.shape[0])
        self.assertEqual(len(rowSumsSel), len(sel_sample))
        i = 0
        for ni in neut_sample:
            num_ones = ni[1].count("1")
            self.assertEqual(num_ones, rowSums[i])
            i += 1

    def test_merge(self):
        m, keys = self.gm.merge()
        s = (
            self.gm_neutral.shape[0] + self.gm_selected.shape[0],
            self.gm_neutral.shape[1],
        )

        self.assertEqual(m.shape, s)

        neutral = np.array(
            [self.pop.mutations[k].neutral for k in keys], dtype=np.int32
        )
        w = np.where(neutral == 1)[0]
        npos = self.gm.neutral.positions
        npos_idx = np.argsort(npos)
        self.assertTrue(np.array_equal(m[w, :], self.gm_neutral[npos_idx, :]))
        self.assertTrue(np.array_equal(self.gm.neutral_matrix[0], m[w, :]))


if __name__ == "__main__":
    unittest.main()
