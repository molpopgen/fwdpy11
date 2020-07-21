import unittest

import numpy as np

import fwdpy11


class testAdditive(unittest.TestCase):
    @classmethod
    def setUp(self):
        GN = fwdpy11.GaussianNoise
        self.w = fwdpy11.Additive(2.0)
        self.t = fwdpy11.Additive(2.0, fwdpy11.GSS(0.0, 1.0))
        self.tn = fwdpy11.Additive(1.0, fwdpy11.GSS(0.0, 1.0), GN(mean=0.1, sd=2.0))

    def testScaling(self):
        self.assertEqual(self.w.scaling, 2.0)
        self.assertEqual(self.t.scaling, 2.0)
        self.assertEqual(self.tn.scaling, 1.0)

    def testFitnessOrTraitBaseClassProperty(self):
        self.assertEqual(self.w.maps_to_fitness, True)
        self.assertEqual(self.w.maps_to_trait_value, False)
        self.assertEqual(self.t.maps_to_fitness, False)
        self.assertEqual(self.t.maps_to_trait_value, True)
        self.assertEqual(self.tn.maps_to_fitness, False)
        self.assertEqual(self.tn.maps_to_trait_value, True)

    def testPickleFitness(self):
        import pickle

        p = pickle.dumps(self.w)
        up = pickle.loads(p)
        self.assertEqual(up.scaling, self.w.scaling)
        self.assertTrue(up.maps_to_fitness)
        self.assertEqual(type(up.noise), type(self.w.noise))
        self.assertEqual(type(up.gvalue_to_fitness), type(self.w.gvalue_to_fitness))

    def testPickleTraitNoNoise(self):
        import pickle

        p = pickle.dumps(self.t)
        up = pickle.loads(p)
        self.assertEqual(up.scaling, self.t.scaling)
        self.assertTrue(up.maps_to_fitness is False)
        self.assertEqual(type(up.noise), type(self.t.noise))
        self.assertEqual(type(up.gvalue_to_fitness), type(self.t.gvalue_to_fitness))

    def testPickleTraitWithNoise(self):
        import pickle

        p = pickle.dumps(self.tn)
        up = pickle.loads(p)
        self.assertEqual(up.scaling, self.tn.scaling)
        self.assertTrue(up.maps_to_fitness is False)
        self.assertEqual(type(up.noise), type(self.tn.noise))
        self.assertEqual(type(up.gvalue_to_fitness), type(self.tn.gvalue_to_fitness))

    def testPickleTraitWithNoiseToFile(self):
        import pickle

        with open("ptest.pickle", "wb") as f:
            pickle.dump(self.tn, f)

        with open("ptest.pickle", "rb") as f:
            up = pickle.load(f)
        self.assertEqual(up.scaling, self.tn.scaling)
        self.assertTrue(up.maps_to_fitness is False)
        self.assertEqual(type(up.noise), type(self.tn.noise))
        self.assertEqual(type(up.gvalue_to_fitness), type(self.tn.gvalue_to_fitness))


class testMultiplicative(unittest.TestCase):
    @classmethod
    def setUp(self):
        GN = fwdpy11.GaussianNoise
        self.w = fwdpy11.Multiplicative(2.0)
        self.t = fwdpy11.Multiplicative(2.0, fwdpy11.GSS(0.0, 1.0))
        self.tn = fwdpy11.Multiplicative(
            1.0, fwdpy11.GSS(0.0, 1.0), GN(mean=0.1, sd=2.0)
        )

    def testScaling(self):
        self.assertEqual(self.w.scaling, 2.0)
        self.assertEqual(self.t.scaling, 2.0)
        self.assertEqual(self.tn.scaling, 1.0)

    def testFitnessOrTraitBaseClassProperty(self):
        self.assertEqual(self.w.maps_to_fitness, True)
        self.assertEqual(self.w.maps_to_trait_value, False)
        self.assertEqual(self.t.maps_to_fitness, False)
        self.assertEqual(self.t.maps_to_trait_value, True)
        self.assertEqual(self.tn.maps_to_fitness, False)
        self.assertEqual(self.tn.maps_to_trait_value, True)

    def testPickleFitness(self):
        import pickle

        p = pickle.dumps(self.w)
        up = pickle.loads(p)
        self.assertEqual(up.scaling, self.w.scaling)
        self.assertTrue(up.maps_to_fitness)
        self.assertEqual(type(up.noise), type(self.w.noise))
        self.assertEqual(type(up.gvalue_to_fitness), type(self.w.gvalue_to_fitness))

    def testPickleTraitNoNoise(self):
        import pickle

        p = pickle.dumps(self.t)
        up = pickle.loads(p)
        self.assertEqual(up.scaling, self.t.scaling)
        self.assertTrue(up.maps_to_fitness is False)
        self.assertEqual(type(up.noise), type(self.t.noise))
        self.assertEqual(type(up.gvalue_to_fitness), type(self.t.gvalue_to_fitness))

    def testPickleTraitWithNoise(self):
        import pickle

        p = pickle.dumps(self.tn)
        up = pickle.loads(p)
        self.assertEqual(up.scaling, self.tn.scaling)
        self.assertTrue(up.maps_to_fitness is False)
        self.assertEqual(type(up.noise), type(self.tn.noise))
        self.assertEqual(type(up.gvalue_to_fitness), type(self.tn.gvalue_to_fitness))


class testGBR(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.gss = fwdpy11.GSS(0.0, 1.0)
        self.gnoise = fwdpy11.GaussianNoise(mean=0.0, sd=1.0)
        self.nonoise = fwdpy11.NoNoise()

    def testPicklingGSS(self):
        import pickle

        gbr = fwdpy11.GBR(self.gss)
        self.assertFalse(gbr.maps_to_fitness)
        self.assertFalse(gbr.maps_to_fitness)
        p = pickle.dumps(gbr, -1)
        up = pickle.loads(p)
        self.assertTrue(up.noise is None)
        self.assertEqual(type(self.gss), type(up.gvalue_to_fitness))

    def testPicklingGSSGaussianNoise(self):
        import pickle

        gbr = fwdpy11.GBR(self.gss, self.gnoise)
        self.assertFalse(gbr.maps_to_fitness)
        self.assertFalse(gbr.maps_to_fitness)
        p = pickle.dumps(gbr, -1)
        up = pickle.loads(p)
        self.assertEqual(type(self.gnoise), type(up.noise))
        self.assertEqual(type(self.gss), type(up.gvalue_to_fitness))


class testGSS(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.x = fwdpy11.GSS(0.0, 1.0)

    def testProperties(self):
        self.assertEqual(self.x.optimum, 0.0)
        self.assertEqual(self.x.VS, 1.0)


if __name__ == "__main__":
    unittest.main()
