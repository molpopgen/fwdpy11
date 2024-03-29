import pickle
import unittest

import fwdpy11
from ezparams import mslike

import custom_additive as ca
import custom_stateless_genotype as general


class testCustomAdditive(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.pop = fwdpy11.DiploidPopulation(1000, 1.0)
        self.pdict = mslike(
            self.pop, dfe=fwdpy11.ExpS(0, 1, 1, -0.05), pneutral=0.95, simlen=10
        )
        self.pdict["prune_selected"] = False
        self.pdict["gvalue"] = ca.additive()
        self.pdict["demography"] = fwdpy11.ForwardDemesGraph.tubes(
            self.pop.deme_sizes()[1], burnin=self.pdict["simlen"], burnin_is_exact=True
        )
        self.rng = fwdpy11.GSLrng(42)
        self.params = fwdpy11.ModelParams(**self.pdict)

    def testEvolve(self):
        fwdpy11.evolvets(self.rng, self.pop, self.params, 10)

    def testPickle(self):
        a = self.params.gvalue
        p = pickle.dumps(a, -1)
        up = pickle.loads(p)
        self.assertEqual(type(a), type(up))


class testGeneralModule(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.pop = fwdpy11.DiploidPopulation(1000, 1.0)
        self.pdict = mslike(
            self.pop,
            dfe=fwdpy11.ConstantS(0, 1, 1, -0.05, 0.05),
            pneutral=0.95,
            simlen=10,
        )
        self.pdict["prune_selected"] = False
        self.pdict["gvalue"] = general.GeneralW()
        self.pdict["demography"] = fwdpy11.ForwardDemesGraph.tubes(
            self.pop.deme_sizes()[1], burnin=self.pdict["simlen"], burnin_is_exact=True
        )
        self.rng = fwdpy11.GSLrng(42)
        self.params = fwdpy11.ModelParams(**self.pdict)

    def testPickle(self):
        a = self.params.gvalue
        p = pickle.dumps(a, -1)
        up = pickle.loads(p)
        self.assertEqual(type(a), type(up))

    def testEvolve(self):
        fwdpy11.evolvets(self.rng, self.pop, self.params, 10)


if __name__ == "__main__":
    unittest.main()
