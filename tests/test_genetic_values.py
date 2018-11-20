import unittest
import fwdpy11
import fwdpy11.genetic_values
import numpy as np


class testSlocusAdditive(unittest.TestCase):
    @classmethod
    def setUp(self):
        from fwdpy11.genetic_value_noise import GaussianNoise
        self.w = fwdpy11.genetic_values.SlocusAdditive(2.0)
        self.t = fwdpy11.genetic_values.SlocusAdditive(
            2.0, fwdpy11.genetic_values.GSS(0.0, 1.0))
        self.tn = fwdpy11.genetic_values.SlocusAdditive(1.0,
                                                        fwdpy11.genetic_values.GSS(
                                                            0.0, 1.0),
                                                        GaussianNoise(mean=0.1, sd=2.0))

    def testScaling(self):
        self.assertEqual(self.w.scaling, 2.0)
        self.assertEqual(self.t.scaling, 2.0)
        self.assertEqual(self.tn.scaling, 1.0)

    def testFitnessOrTrait(self):
        self.assertEqual(self.w.is_fitness, True)
        self.assertEqual(self.t.is_fitness, False)
        self.assertEqual(self.tn.is_fitness, False)

    def testPickleFitness(self):
        import pickle
        p = pickle.dumps(self.w)
        up = pickle.loads(p)
        self.assertEqual(up.scaling, self.w.scaling)
        self.assertTrue(up.is_fitness)
        self.assertEqual(type(up.noise), type(self.w.noise))
        self.assertEqual(type(up.gvalue_to_fitness),
                         type(self.w.gvalue_to_fitness))

    def testPickleTraitNoNoise(self):
        import pickle
        p = pickle.dumps(self.t)
        up = pickle.loads(p)
        self.assertEqual(up.scaling, self.t.scaling)
        self.assertTrue(up.is_fitness is False)
        self.assertEqual(type(up.noise), type(self.t.noise))
        self.assertEqual(type(up.gvalue_to_fitness),
                         type(self.t.gvalue_to_fitness))

    def testPickleTraitWithNoise(self):
        import pickle
        p = pickle.dumps(self.tn)
        up = pickle.loads(p)
        self.assertEqual(up.scaling, self.tn.scaling)
        self.assertTrue(up.is_fitness is False)
        self.assertEqual(type(up.noise), type(self.tn.noise))
        self.assertEqual(type(up.gvalue_to_fitness),
                         type(self.tn.gvalue_to_fitness))

    def testPickleTraitWithNoiseToFile(self):
        import pickle
        with open("ptest.pickle", "wb") as f:
            pickle.dump(self.tn, f)

        with open("ptest.pickle", "rb") as f:
            up = pickle.load(f)
        self.assertEqual(up.scaling, self.tn.scaling)
        self.assertTrue(up.is_fitness is False)
        self.assertEqual(type(up.noise), type(self.tn.noise))
        self.assertEqual(type(up.gvalue_to_fitness),
                         type(self.tn.gvalue_to_fitness))


class testMlocusAdditive(unittest.TestCase):
    @classmethod
    def setUp(self):
        from fwdpy11.genetic_value_noise import GaussianNoise
        self.w = fwdpy11.genetic_values.MlocusAdditive(2.0)
        self.t = fwdpy11.genetic_values.MlocusAdditive(
            2.0, fwdpy11.genetic_values.GSS(0.0, 1.0))
        self.tn = fwdpy11.genetic_values.MlocusAdditive(1.0,
                                                        fwdpy11.genetic_values.GSS(
                                                            0.0, 1.0),
                                                        GaussianNoise(mean=0.1, sd=2.0))

    def testScaling(self):
        self.assertEqual(self.w.scaling, 2.0)
        self.assertEqual(self.t.scaling, 2.0)
        self.assertEqual(self.tn.scaling, 1.0)

    def testFitnessOrTrait(self):
        self.assertEqual(self.w.is_fitness, True)
        self.assertEqual(self.t.is_fitness, False)
        self.assertEqual(self.tn.is_fitness, False)

    def testPickleFitness(self):
        import pickle
        p = pickle.dumps(self.w)
        up = pickle.loads(p)
        self.assertEqual(up.scaling, self.w.scaling)
        self.assertTrue(up.is_fitness)
        self.assertEqual(type(up.noise), type(self.w.noise))
        self.assertEqual(type(up.gvalue_to_fitness),
                         type(self.w.gvalue_to_fitness))

    def testPickleTraitNoNoise(self):
        import pickle
        p = pickle.dumps(self.t)
        up = pickle.loads(p)
        self.assertEqual(up.scaling, self.t.scaling)
        self.assertTrue(up.is_fitness is False)
        self.assertEqual(type(up.noise), type(self.t.noise))
        self.assertEqual(type(up.gvalue_to_fitness),
                         type(self.t.gvalue_to_fitness))

    def testPickleTraitWithNoise(self):
        import pickle
        p = pickle.dumps(self.tn)
        up = pickle.loads(p)
        self.assertEqual(up.scaling, self.tn.scaling)
        self.assertTrue(up.is_fitness is False)
        self.assertEqual(type(up.noise), type(self.tn.noise))
        self.assertEqual(type(up.gvalue_to_fitness),
                         type(self.tn.gvalue_to_fitness))


class testSlocusMult(unittest.TestCase):
    @classmethod
    def setUp(self):
        from fwdpy11.genetic_value_noise import GaussianNoise
        self.w = fwdpy11.genetic_values.SlocusMult(2.0)
        self.t = fwdpy11.genetic_values.SlocusMult(
            2.0, fwdpy11.genetic_values.GSS(0.0, 1.0))
        self.tn = fwdpy11.genetic_values.SlocusMult(1.0,
                                                    fwdpy11.genetic_values.GSS(
                                                        0.0, 1.0),
                                                    GaussianNoise(mean=0.1, sd=2.0))

    def testScaling(self):
        self.assertEqual(self.w.scaling, 2.0)
        self.assertEqual(self.t.scaling, 2.0)
        self.assertEqual(self.tn.scaling, 1.0)

    def testFitnessOrTrait(self):
        self.assertEqual(self.w.is_fitness, True)
        self.assertEqual(self.t.is_fitness, False)
        self.assertEqual(self.tn.is_fitness, False)

    def testPickleFitness(self):
        import pickle
        p = pickle.dumps(self.w)
        up = pickle.loads(p)
        self.assertEqual(up.scaling, self.w.scaling)
        self.assertTrue(up.is_fitness)
        self.assertEqual(type(up.noise), type(self.w.noise))
        self.assertEqual(type(up.gvalue_to_fitness),
                         type(self.w.gvalue_to_fitness))

    def testPickleTraitNoNoise(self):
        import pickle
        p = pickle.dumps(self.t)
        up = pickle.loads(p)
        self.assertEqual(up.scaling, self.t.scaling)
        self.assertTrue(up.is_fitness is False)
        self.assertEqual(type(up.noise), type(self.t.noise))
        self.assertEqual(type(up.gvalue_to_fitness),
                         type(self.t.gvalue_to_fitness))

    def testPickleTraitWithNoise(self):
        import pickle
        p = pickle.dumps(self.tn)
        up = pickle.loads(p)
        self.assertEqual(up.scaling, self.tn.scaling)
        self.assertTrue(up.is_fitness is False)
        self.assertEqual(type(up.noise), type(self.tn.noise))
        self.assertEqual(type(up.gvalue_to_fitness),
                         type(self.tn.gvalue_to_fitness))


class testMlocusMult(unittest.TestCase):
    @classmethod
    def setUp(self):
        from fwdpy11.genetic_value_noise import GaussianNoise
        self.w = fwdpy11.genetic_values.MlocusMult(2.0)
        self.t = fwdpy11.genetic_values.MlocusMult(
            2.0, fwdpy11.genetic_values.GSS(0.0, 1.0))
        self.tn = fwdpy11.genetic_values.MlocusMult(1.0,
                                                    fwdpy11.genetic_values.GSS(
                                                        0.0, 1.0),
                                                    GaussianNoise(mean=0.1, sd=2.0))

    def testScaling(self):
        self.assertEqual(self.w.scaling, 2.0)
        self.assertEqual(self.t.scaling, 2.0)
        self.assertEqual(self.tn.scaling, 1.0)

    def testFitnessOrTrait(self):
        self.assertEqual(self.w.is_fitness, True)
        self.assertEqual(self.t.is_fitness, False)
        self.assertEqual(self.tn.is_fitness, False)

    def testPickleFitness(self):
        import pickle
        p = pickle.dumps(self.w)
        up = pickle.loads(p)
        self.assertEqual(up.scaling, self.w.scaling)
        self.assertTrue(up.is_fitness)
        self.assertEqual(type(up.noise), type(self.w.noise))
        self.assertEqual(type(up.gvalue_to_fitness),
                         type(self.w.gvalue_to_fitness))

    def testPickleTraitNoNoise(self):
        import pickle
        p = pickle.dumps(self.t)
        up = pickle.loads(p)
        self.assertEqual(up.scaling, self.t.scaling)
        self.assertTrue(up.is_fitness is False)
        self.assertEqual(type(up.noise), type(self.t.noise))
        self.assertEqual(type(up.gvalue_to_fitness),
                         type(self.t.gvalue_to_fitness))

    def testPickleTraitWithNoise(self):
        import pickle
        p = pickle.dumps(self.tn)
        up = pickle.loads(p)
        self.assertEqual(up.scaling, self.tn.scaling)
        self.assertTrue(up.is_fitness is False)
        self.assertEqual(type(up.noise), type(self.tn.noise))
        self.assertEqual(type(up.gvalue_to_fitness),
                         type(self.tn.gvalue_to_fitness))


class testSlocusGBR(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.gss = fwdpy11.genetic_values.GSS(0.0, 1.0)
        self.gnoise = fwdpy11.genetic_value_noise.GaussianNoise(
            mean=0.0, sd=1.0)
        self.nonoise = fwdpy11.genetic_value_noise.NoNoise()

    def testPicklingGSS(self):
        import pickle
        gbr = fwdpy11.genetic_values.SlocusGBR(self.gss)
        p = pickle.dumps(gbr, -1)
        up = pickle.loads(p)
        self.assertEqual(type(self.nonoise), type(up.noise))
        self.assertEqual(type(self.gss), type(up.gvalue_to_fitness))

    def testPicklingGSSGaussianNoise(self):
        import pickle
        gbr = fwdpy11.genetic_values.SlocusGBR(self.gss, self.gnoise)
        p = pickle.dumps(gbr, -1)
        up = pickle.loads(p)
        self.assertEqual(type(self.gnoise), type(up.noise))
        self.assertEqual(type(self.gss), type(up.gvalue_to_fitness))


class testMlocusGBR(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.gss = fwdpy11.genetic_values.GSS(0.0, 1.0)
        self.gnoise = fwdpy11.genetic_value_noise.GaussianNoise(
            mean=0.0, sd=1.0)
        self.nonoise = fwdpy11.genetic_value_noise.NoNoise()

    def testPicklingGSS(self):
        import pickle
        gbr = fwdpy11.genetic_values.MlocusGBR(self.gss)
        p = pickle.dumps(gbr, -1)
        up = pickle.loads(p)
        self.assertEqual(type(self.nonoise), type(up.noise))
        self.assertEqual(type(self.gss), type(up.gvalue_to_fitness))

    def testPicklingGSSGaussianNoise(self):
        import pickle
        gbr = fwdpy11.genetic_values.MlocusGBR(self.gss, self.gnoise)
        p = pickle.dumps(gbr, -1)
        up = pickle.loads(p)
        self.assertEqual(type(self.gnoise), type(up.noise))
        self.assertEqual(type(self.gss), type(up.gvalue_to_fitness))


class testGSS(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.x = fwdpy11.genetic_values.GSS(0.0, 1.0)

    def testProperties(self):
        self.assertEqual(self.x.opt, 0.0)
        self.assertEqual(self.x.VS, 1.0)


class testGSSandGSSmoConsistency(unittest.TestCase):
    """
    This tests that GSS and GSSmo have opt and VS
    in the same order.  If that were not true,
    fitness calculations would come out differently
    and the test would fail.
    """
    @classmethod
    def setUp(self):
        self.a = fwdpy11.genetic_values.SlocusAdditive(
            2.0, fwdpy11.genetic_values.GSS(0.0, 1.0))
        self.b = fwdpy11.genetic_values.SlocusAdditive(
            2.0, fwdpy11.genetic_values.GSSmo([(0, 0.0, 1.0)]))
        self.pop = fwdpy11.SlocusPop(1000)

    def testFitnesses(self):
        wa = [self.a.fitness(i, self.pop) for i in range(self.pop.N)]
        wb = [self.b.fitness(i, self.pop) for i in range(self.pop.N)]
        for i, j in zip(wa, wb):
            self.assertEqual(i, j)


if __name__ == "__main__":
    unittest.main()
