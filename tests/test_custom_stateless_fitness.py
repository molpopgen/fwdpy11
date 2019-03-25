import fwdpy11
import fwdpy11.ezparams
import fwdpy11
import pickle
import unittest
import custom_additive as ca
import custom_stateless_genotype as general


class testCustomAdditive(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.pop = fwdpy11.DiploidPopulation(1000)
        self.pdict = fwdpy11.ezparams.mslike(self.pop,
                                             dfe=fwdpy11.ExpS(0, 1, 1, -0.05),
                                             pneutral=0.95, simlen=10)
        self.pdict['gvalue'] = ca.additive()
        self.rng = fwdpy11.GSLrng(42)
        self.params = fwdpy11.ModelParams(**self.pdict)

    def testEvolve(self):
        fwdpy11.evolve_genomes(self.rng, self.pop, self.params)

    def testPickle(self):
        a = self.params.gvalue
        p = pickle.dumps(a, -1)
        up = pickle.loads(p)
        self.assertEqual(type(a), type(up))

    # TODO: test this once built-in DiploidAdditive is callable
    # def testCorrectNess(self):
    #     fwdpy11.evolve_genomes(self.rng, self.pop, self.params)
    #     a = fwdpy11.fitness.DiploidAdditive(2.0)
    #     for i in self.pop.diploids:
    #         self.assertEqual(i.w, a(i, self.pop))


class testGeneralModule(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.pop = fwdpy11.DiploidPopulation(1000)
        self.pdict = fwdpy11.ezparams.mslike(self.pop,
                                             dfe=fwdpy11.ConstantS(
                                                 0, 1, 1, -0.05, 0.05),
                                             pneutral=0.95, simlen=10)
        self.pdict['gvalue'] = general.GeneralW()
        self.rng = fwdpy11.GSLrng(42)
        self.params = fwdpy11.ModelParams(**self.pdict)

    def testPickle(self):
        a = self.params.gvalue
        p = pickle.dumps(a, -1)
        up = pickle.loads(p)
        self.assertEqual(type(a), type(up))

    def testEvolve(self):
        fwdpy11.evolve_genomes(self.rng, self.pop, self.params)


if __name__ == "__main__":
    unittest.main()
