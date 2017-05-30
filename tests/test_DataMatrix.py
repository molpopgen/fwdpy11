import unittest
import fwdpy11 as fp11
import fwdpy11.sampling
import numpy as np
from quick_pops import quick_nonneutral_slocus

class test_DataMatrixFromSlocusPop(unittest.TestCase):
    """
    Much of this is already tested in fwdpp's unit tests,
    but it never hurts to try another RNG seed, etc.
    """
    @classmethod
    def setUpClass(self):
        self.pop = quick_nonneutral_slocus()
        self.indlist = [i for i in range(100,150)]
        self.keys = fwdpy11.sampling.mutation_keys(self.pop,self.indlist)
        self.hm = fwdpy11.sampling.haplotype_matrix(self.pop,self.indlist,self.keys[0],self.keys[1])
        self.gm = fwdpy11.sampling.genotype_matrix(self.pop,self.indlist,self.keys[0],self.keys[1])
    def testKeyNeutralityAndCount(self):
        for i in self.keys[0]:
            self.assertTrue(self.pop.mutations[i[0]].neutral)
            self.assertTrue(self.pop.mcounts[i[0]]>0)
        for i in self.keys[1]:
            self.assertFalse(self.pop.mutations[i[0]].neutral)
            self.assertTrue(self.pop.mcounts[i[0]]>0)
    def testHapMat(self):
        nm = self.hm.neutral()
        self.assertTrue(nm.dtype == np.int8)
        #convert from 1d to 2d (future versions will simply return a 
        #2d matrix, but we need a new pybind11 release for that):
        nm=nm.reshape((self.hm.nrow,int(len(nm)/self.hm.nrow)))
        self.assertEqual(nm.ndim,2)
        sm = self.hm.selected()
        sm=sm.reshape((self.hm.nrow,int(len(sm)/self.hm.nrow)))
        self.assertTrue(sm.dtype == np.int8)
        #Get the row sums
        rowSums = nm.sum(axis=1)
        rowSumsSel = sm.sum(axis=1)
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
        nm = self.hm.neutral()
        self.assertTrue(nm.dtype == np.int8)
        #convert from 1d to 2d (future versions will simply return a 
        #2d matrix, but we need a new pybind11 release for that):
        nm=nm.reshape((self.gm.nrow,int(len(nm)/self.gm.nrow)))
        self.assertEqual(nm.ndim,2)
        sm = self.gm.selected()
        sm=sm.reshape((self.gm.nrow,int(len(sm)/self.gm.nrow)))
        self.assertTrue(sm.dtype == np.int8)
        #Get the row sums
        rowSums = nm.sum(axis=1)
        rowSumsSel = sm.sum(axis=1)
        self.assertEqual(len(rowSums),self.gm.nrow)
        self.assertEqual(len(rowSumsSel),self.gm.nrow)
        j=0
        for i in range(100,150):
            #Num neutral variants in diploid i, gamete 0
            nmuts = len(self.pop.gametes[self.pop.diploids[i].first].mutations)
            nmuts += len(self.pop.gametes[self.pop.diploids[i].second].mutations)
            self.assertEqual(nmuts,rowSums[j])

            #Now, test numbers of selected
            nmuts = len(self.pop.gametes[self.pop.diploids[i].first].smutations)
            nmuts += len(self.pop.gametes[self.pop.diploids[i].second].smutations)
            self.assertEqual(nmuts,rowSumsSel[j])

            j+=1
