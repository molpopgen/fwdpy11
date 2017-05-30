import unittest
import fwdpy11 as fp11
import fwdpy11.sampling
import numpy as np
from quick_pops import quick_nonneutral_slocus

class test_DataMatrixFromSlocusPop(unittest.TestCase):
    """
    Much of this is already tested in fwdpp's unit tests,
    but the purpose here is to make sure that the wrappers
    are working appropriately.
    """
    @classmethod
    def setUpClass(self):
        self.pop = quick_nonneutral_slocus()
        self.indlist = [i for i in range(100,150)]
        self.keys = fwdpy11.sampling.mutation_keys(self.pop,self.indlist)
        self.hm = fwdpy11.sampling.haplotype_matrix(self.pop,self.indlist,self.keys[0],self.keys[1])
        self.gm = fwdpy11.sampling.genotype_matrix(self.pop,self.indlist,self.keys[0],self.keys[1])
        self.hm_neutral = self.hm.neutral()
        self.hm_neutral = self.hm_neutral.reshape((self.hm.nrow,int(len(self.hm_neutral)/self.hm.nrow)))
        self.hm_selected = self.hm.selected()
        self.hm_selected = self.hm_selected.reshape((self.hm.nrow,int(len(self.hm_selected)/self.hm.nrow)))
        self.gm_neutral = self.gm.neutral()
        self.gm_neutral = self.gm_neutral.reshape((self.gm.nrow,int(len(self.gm_neutral)/self.gm.nrow)))
        self.gm_selected = self.gm.selected()
        self.gm_selected = self.gm_selected.reshape((self.gm.nrow,int(len(self.gm_selected)/self.gm.nrow)))
    def testKeyNeutralityAndCount(self):
        for i in self.keys[0]:
            self.assertTrue(self.pop.mutations[i[0]].neutral)
            self.assertTrue(self.pop.mcounts[i[0]]>0)
        for i in self.keys[1]:
            self.assertFalse(self.pop.mutations[i[0]].neutral)
            self.assertTrue(self.pop.mcounts[i[0]]>0)
    def testHapMat(self):
        self.assertTrue(self.hm_neutral.dtype == np.int8)
        self.assertEqual(self.hm_neutral.ndim,2)
        self.assertTrue(self.hm_selected.dtype == np.int8)
        self.assertEqual(self.hm_selected.ndim,2)
        #Get the row sums
        rowSums = self.hm_neutral.sum(axis=1)
        rowSumsSel = self.hm_selected.sum(axis=1)
        self.assertEqual(len(rowSums),self.hm.nrow)
        self.assertEqual(len(rowSumsSel),self.hm.nrow)
        j=0
        for i in range(100,150):
            #Num neutral variants in diploid i, gamete 0
            nmuts = len(self.pop.gametes[self.pop.diploids[i].first].mutations)
            self.assertEqual(nmuts,rowSums[j])
            #Num neutral variants in diploid i, gamete 1
            nmuts = len(self.pop.gametes[self.pop.diploids[i].second].mutations)
            self.assertEqual(nmuts,rowSums[j+1])

            #Now, test numbers of selected
            nmuts = len(self.pop.gametes[self.pop.diploids[i].first].smutations)
            self.assertEqual(nmuts,rowSumsSel[j])
            nmuts = len(self.pop.gametes[self.pop.diploids[i].second].smutations)
            self.assertEqual(nmuts,rowSumsSel[j+1])

            j+=2
    def testGenoMat(self):
        self.assertTrue(self.gm_neutral.dtype == np.int8)
        self.assertEqual(self.gm_neutral.ndim,2)
        #Get the row sums
        rowSums = self.gm_neutral.sum(axis=1)
        rowSumsSel = self.gm_selected.sum(axis=1)
        self.assertEqual(len(rowSums),self.gm.nrow)
        self.assertEqual(len(rowSumsSel),self.gm.nrow)
        j=0
        for i in range(100,150):
            nmuts = len(self.pop.gametes[self.pop.diploids[i].first].mutations)
            nmuts += len(self.pop.gametes[self.pop.diploids[i].second].mutations)
            self.assertEqual(nmuts,rowSums[j])

            #Now, test numbers of selected
            nmuts = len(self.pop.gametes[self.pop.diploids[i].first].smutations)
            nmuts += len(self.pop.gametes[self.pop.diploids[i].second].smutations)
            self.assertEqual(nmuts,rowSumsSel[j])

            j+=1
    def testHapMatToSample(self):
        neut_sample = fwdpy11.sampling.matrix_to_sample(self.hm,True)
        sel_sample = fwdpy11.sampling.matrix_to_sample(self.hm,False)
        colSums = self.hm_neutral.sum(axis=0)
        colSumsSel = self.hm_selected.sum(axis=0)
        self.assertEqual(len(neut_sample),self.hm_neutral.shape[1])
        self.assertEqual(len(colSums),len(neut_sample))
        self.assertEqual(len(sel_sample),self.hm_selected.shape[1])
        self.assertEqual(len(colSumsSel),len(sel_sample))
        i=0
        for ni in neut_sample:
            num_ones = ni[1].count('1')
            self.assertEqual(num_ones,colSums[i])
            i+=1
