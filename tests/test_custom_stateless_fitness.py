import cppimport
import fwdpy11
import fwdpy11.ezparams
import fwdpy11.model_params
import fwdpy11.wright_fisher
import fwdpy11.fwdpy11_types
import pickle
import unittest
cppimport.force_rebuild()
ca = cppimport.imp("custom_additive")
general = cppimport.imp("custom_stateless_genotype")


class testCustomAdditive(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.pop = fwdpy11.SlocusPop(1000)
        self.pdict = fwdpy11.ezparams.mslike(self.pop,
                                             dfe=fwdpy11.ExpS(0, 1, 1, -0.05),
                                             pneutral=0.95, simlen=10)
        self.pdict['gvalue'] = ca.additive()
        self.rng = fwdpy11.GSLrng(42)
        self.params = fwdpy11.model_params.ModelParams(**self.pdict)

    def testEvolve(self):
        fwdpy11.wright_fisher.evolve(self.rng, self.pop, self.params)

    def testPickle(self):
        a = self.params.gvalue
        p = pickle.dumps(a, -1)
        up = pickle.loads(p)
        self.assertEqual(type(a), type(up))

    # TODO: test this once built-in SlocusAdditive is callable
    # def testCorrectNess(self):
    #     fwdpy11.wright_fisher.evolve(self.rng, self.pop, self.params)
    #     a = fwdpy11.fitness.SlocusAdditive(2.0)
    #     for i in self.pop.diploids:
    #         self.assertEqual(i.w, a(i, self.pop))


class testGeneralModule(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.pop = fwdpy11.SlocusPop(1000)
        self.pdict = fwdpy11.ezparams.mslike(self.pop,
                                             dfe=fwdpy11.ConstantS(
                                                 0, 1, 1, -0.05, 0.05),
                                             pneutral=0.95, simlen=10)
        self.pdict['gvalue'] = general.GeneralW()
        self.rng = fwdpy11.GSLrng(42)
        self.params = fwdpy11.model_params.ModelParams(**self.pdict)

    def testPickle(self):
        a = self.params.gvalue
        p = pickle.dumps(a, -1)
        up = pickle.loads(p)
        self.assertEqual(type(a), type(up))

    def testEvolve(self):
        fwdpy11.wright_fisher.evolve(self.rng, self.pop, self.params)


if __name__ == "__main__":
    unittest.main()
