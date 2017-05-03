import unittest
import fwdpy11 as fp11
import fwdpy11.wright_fisher as wf
import numpy as np

class testWFevolve(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.pop = fp11.SlocusPop(1000)
        self.rng = fp11.GSLrng(42)
    def testEvolve(self):
        p = fp11.model_params.SlocusParams()
        p.rates=(1e-3,1e-3,1e-3)
        p.demography=np.array([1000]*100,dtype=np.uint32)
        p.nregions=[fp11.Region(0,1,1)]
        p.sregions=[fp11.ExpS(0,1,1,-1e-2)]
        p.recregions=p.nregions
        from fwdpy11.wright_fisher import evolve
        evolve(self.rng,self.pop,p)

