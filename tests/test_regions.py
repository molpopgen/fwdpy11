# Unit tests for fwdpy11.regions

import pickle
import unittest

import numpy as np

import call_Sregion
import fwdpy11

# Choose global values that differ from default values in Python class constructors
BEG, END, WEIGHT, DOM, LABEL, COUPLED = 0.0, 1.0, 0.5, 0.5, 63, False


class testRegion(unittest.TestCase):
    """
    These tests all apply to Sregions b/c an
    Sregion contains a Region
    """

    def testBadBoundaries(self):
        with self.assertRaises(ValueError):
            fwdpy11.Region(np.nan, 1, 1.0)
        with self.assertRaises(ValueError):
            fwdpy11.Region(0, np.nan, 1.0)
        with self.assertRaises(ValueError):
            fwdpy11.Region(0, 0, 1.0)
        with self.assertRaises(ValueError):
            fwdpy11.Region(0, -1, 1.0)

    def testBadWeight(self):
        with self.assertRaises(ValueError):
            fwdpy11.Region(0, 1, np.nan)
        with self.assertRaises(ValueError):
            fwdpy11.Region(0, 1, -1e-3)


class test_PickleRegion(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.r = fwdpy11.Region(BEG, END, WEIGHT, COUPLED, LABEL)

    def test_pickling(self):
        p = pickle.dumps(self.r, -1)
        up = pickle.loads(p)
        self.assertEqual(self.r.beg, up.beg)
        self.assertEqual(self.r.end, up.end)
        self.assertEqual(self.r.weight, up.weight)
        self.assertEqual(self.r.coupled, up.coupled)
        self.assertEqual(self.r.label, up.label)


class testExpS(unittest.TestCase):
    # NOTE: this tests applies to all Sregions
    # b/c the scaling param is part of the base class
    def testBadScaling(self):
        with self.assertRaises(ValueError):
            fwdpy11.ExpS(0, 1, 1.0, -0.1, scaling=np.nan)

    def testBadMean(self):
        with self.assertRaises(ValueError):
            fwdpy11.ExpS(0, 1, 1, np.nan)

    def testBadDominance(self):
        with self.assertRaises(ValueError):
            fwdpy11.ExpS(0, 1, 1, -1.0, h=np.nan)


class test_PickleExpS(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.r = fwdpy11.ExpS(BEG, END, WEIGHT, -0.2, DOM, COUPLED, LABEL)

    def test_pickling(self):
        p = pickle.dumps(self.r, -1)
        up = pickle.loads(p)
        self.assertEqual(self.r.beg, up.beg)
        self.assertEqual(self.r.end, up.end)
        self.assertEqual(self.r.weight, up.weight)
        self.assertEqual(self.r.coupled, up.coupled)
        self.assertEqual(self.r.label, up.label)
        self.assertEqual(up.mean, -0.2)
        self.assertEqual(up.h, DOM)


class testGammaS(unittest.TestCase):
    def testBadMean(self):
        with self.assertRaises(ValueError):
            fwdpy11.GammaS(0, 1, 1, mean=np.nan, shape_parameter=1.0)

    def testBadShape(self):
        with self.assertRaises(ValueError):
            fwdpy11.GammaS(0, 1, 1, mean=1.0, shape_parameter=np.nan)

    def testBadDominance(self):
        with self.assertRaises(ValueError):
            fwdpy11.GammaS(0, 1, 1, mean=1.0, shape_parameter=1.0, h=np.nan)


class test_PickleGammaS(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.r = fwdpy11.GammaS(BEG, END, WEIGHT, -0.2, 0.5, DOM, COUPLED, LABEL)

    def test_pickling(self):
        p = pickle.dumps(self.r, -1)
        up = pickle.loads(p)
        self.assertEqual(self.r.beg, up.beg)
        self.assertEqual(self.r.end, up.end)
        self.assertEqual(self.r.weight, up.weight)
        self.assertEqual(self.r.coupled, up.coupled)
        self.assertEqual(self.r.label, up.label)
        self.assertEqual(up.mean, -0.2)
        self.assertEqual(up.shape_parameter, 0.5)
        self.assertEqual(up.h, DOM)


class testUniformS(unittest.TestCase):
    def testBadRange(self):
        with self.assertRaises(ValueError):
            fwdpy11.UniformS(0, 1, 1, 0, 0)
        with self.assertRaises(ValueError):
            fwdpy11.UniformS(0, 1, 1, np.nan, 1)
        with self.assertRaises(ValueError):
            fwdpy11.UniformS(0, 1, 1, 1, np.nan)

    def testBadDominance(self):
        with self.assertRaises(ValueError):
            fwdpy11.UniformS(0, 1, 1, lo=1.0, hi=1.0, h=np.nan)


class test_PickleUniformS(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.r = fwdpy11.UniformS(BEG, END, WEIGHT, -0.2, 0.5, DOM, COUPLED, LABEL)

    def test_pickling(self):
        p = pickle.dumps(self.r, -1)
        up = pickle.loads(p)
        self.assertEqual(self.r.beg, up.beg)
        self.assertEqual(self.r.end, up.end)
        self.assertEqual(self.r.weight, up.weight)
        self.assertEqual(self.r.coupled, up.coupled)
        self.assertEqual(self.r.label, up.label)
        self.assertEqual(up.lo, -0.2)
        self.assertEqual(up.hi, 0.5)
        self.assertEqual(up.h, DOM)


class testGaussianS(unittest.TestCase):
    def testBadStdDev(self):
        with self.assertRaises(ValueError):
            fwdpy11.GaussianS(0, 1, 1, np.nan)
        with self.assertRaises(ValueError):
            fwdpy11.GaussianS(0, 1, 1, 0)
        with self.assertRaises(ValueError):
            fwdpy11.GaussianS(0, 1, 1, -1e-3)

    def testBadDominance(self):
        with self.assertRaises(ValueError):
            fwdpy11.GaussianS(0, 1, 1, 1e-3, h=np.nan)


class test_PickleGaussianS(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.r = fwdpy11.GaussianS(BEG, END, WEIGHT, 0.5, DOM, COUPLED, LABEL)

    def test_pickling(self):
        p = pickle.dumps(self.r, -1)
        up = pickle.loads(p)
        self.assertEqual(self.r.beg, up.beg)
        self.assertEqual(self.r.end, up.end)
        self.assertEqual(self.r.weight, up.weight)
        self.assertEqual(self.r.coupled, up.coupled)
        self.assertEqual(self.r.label, up.label)
        self.assertEqual(up.sd, 0.5)
        self.assertEqual(up.h, DOM)


class testConstantS(unittest.TestCase):
    def testBadValue(self):
        with self.assertRaises(ValueError):
            fwdpy11.ConstantS(0, 1, 1, np.nan)

    def testBadDominance(self):
        with self.assertRaises(ValueError):
            fwdpy11.ConstantS(0, 1, 1, 1, np.nan)


class test_PickleConstantS(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.r = fwdpy11.ConstantS(BEG, END, WEIGHT, 0.5, DOM, COUPLED, LABEL)

    def test_pickling(self):
        p = pickle.dumps(self.r, -1)
        up = pickle.loads(p)
        self.assertEqual(self.r.beg, up.beg)
        self.assertEqual(self.r.end, up.end)
        self.assertEqual(self.r.weight, up.weight)
        self.assertEqual(self.r.coupled, up.coupled)
        self.assertEqual(self.r.label, up.label)
        self.assertEqual(up.s, 0.5)
        self.assertEqual(up.h, DOM)


class Test_mvDES(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.r = fwdpy11.mvDES(
            [fwdpy11.ExpS(0, 1, 1, 0.1), fwdpy11.ExpS(0, 1, 1, -0.1)],
            np.zeros(2),
            np.identity(2),
        )

    def test_pickling(self):
        p = pickle.dumps(self.r, -1)
        up = pickle.loads(p)
        self.assertTrue(np.array_equal(up.means, self.r.means))
        self.assertTrue(np.array_equal(up.matrix, np.identity(2)))
        self.assertTrue(np.array_equal(up.dominance, self.r.dominance))

    def test_bad_init_1(self):
        with self.assertRaises(ValueError):
            fwdpy11.mvDES(
                [fwdpy11.ExpS(0, 1, 1, 0.1), fwdpy11.ExpS(0, 1, 1, -0.1)],
                np.zeros(3),
                np.identity(2),
            )

    def test_mvLogNormalS(self):
        try:
            mvln = fwdpy11.LogNormalS.mv(0, 1, 1)
            fwdpy11.mvDES(mvln, np.zeros(5), np.identity(5))
        except:  # NOQA
            self.fail("unexpected exception")

    def test_mvLogNormalS_pickling(self):
        mvln = fwdpy11.LogNormalS.mv(0, 1, 1)
        des = fwdpy11.mvDES(mvln, np.zeros(5), np.identity(5))
        p = pickle.dumps(des, -1)
        pickle.loads(p)

    def test_mvLogNormalS_bad_init(self):
        mvln = fwdpy11.LogNormalS.mv(0, 1, 1)
        with self.assertRaises(ValueError):
            fwdpy11.mvDES([mvln], np.zeros(5), np.identity(5))

    def test_MultivariateGaussian(self):
        try:
            mvg = fwdpy11.MultivariateGaussianEffects(0, 1, 1, np.identity(4))
            fwdpy11.mvDES(mvg, np.zeros(4))
        except:  # NOQA
            self.fail("unexpected exception")

    def test_MultivariateGaussian_pickling(self):
        mvg = fwdpy11.MultivariateGaussianEffects(0, 1, 1, np.identity(4))
        des = fwdpy11.mvDES(mvg, np.zeros(4))
        p = pickle.dumps(des, -1)
        pickle.loads(p)

    def test_MultivariateGaussian_bad_init(self):
        mvg = fwdpy11.MultivariateGaussianEffects(0, 1, 1, np.identity(4))
        with self.assertRaises(ValueError):
            fwdpy11.mvDES([mvg], np.zeros(4), np.identity(4))


class testMutationRegions(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.nregions = [fwdpy11.Region(0, 1, 1)]
        self.sregions = [
            fwdpy11.ConstantS(0, 1, 1, -0.1, 0.25),
            fwdpy11.ExpS(0, 1, 0.001, 0.01, 1.0),
        ]

    def test_create(self):
        mr = fwdpy11.MutationRegions.create(0.5, self.nregions, self.sregions)  # NOQA

    def test_shape(self):
        self.assertEqual(self.sregions[0].shape, (1,))


class testMultivariateGaussianEffects(unittest.TestCase):
    def testConstruction(self):
        try:
            fwdpy11.MultivariateGaussianEffects(0, 1, 1, np.identity(2))
        except Exception:
            self.fail("Unexpected exception during object construction")

    def test_shape(self):
        mv = fwdpy11.MultivariateGaussianEffects(0, 1, 1, np.identity(2))
        self.assertEqual(mv.shape, (2,))

    def testInvalidMatrixShape(self):
        # Input matrix has wrong shape:
        m = np.array([-1.0] * 4)
        with self.assertRaises(ValueError):
            fwdpy11.MultivariateGaussianEffects(0, 1, 1, m)

    def testInvalidMatrixShape2(self):
        m = np.array([-1.0] * 6).reshape((3, 2))
        with self.assertRaises(ValueError):
            fwdpy11.MultivariateGaussianEffects(0, 1, 1, m)

    def testInvalidMatrix(self):
        # This is not a positive-definite covariance
        # matrix:
        m = np.array([-1.0] * 4).reshape((2, 2))
        with self.assertRaises(ValueError):
            fwdpy11.MultivariateGaussianEffects(0, 1, 1, m)

    def testMatrixWithNAN(self):
        m = np.identity(4).reshape((4, 4))
        m[0, 1] = np.nan
        with self.assertRaises(ValueError):
            fwdpy11.MultivariateGaussianEffects(0, 1, 1, m)


class test_PickleMultivariateGaussianEffects(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.r = fwdpy11.MultivariateGaussianEffects(
            BEG, END, WEIGHT, np.identity(2), -0.3, DOM, COUPLED, LABEL
        )

    def test_pickling(self):
        p = pickle.dumps(self.r, -1)
        up = pickle.loads(p)
        self.assertEqual(self.r, up)


class testRecombinationRegions(unittest.TestCase):
    def test_create(self):
        rr = fwdpy11.RecombinationRegions(
            1e-3, [fwdpy11.Region(0, 1, 1), fwdpy11.Region(1, 2, 1)]
        )
        self.assertEqual(len(rr.weights), 2)

    def test_bad_input(self):
        with self.assertRaises(TypeError):
            fwdpy11.RecombinationRegions(
                1e-3, [fwdpy11.ExpS(0, 1, 1, -0.2), fwdpy11.GaussianS(1, 2, 1, 0.25)]
            )


class TestLogNormalS(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.lns = fwdpy11.LogNormalS(0, 1, 1, 0.1, 2)

    def test_create(self):
        self.assertEqual(self.lns.beg, 0)
        self.assertEqual(self.lns.end, 1)
        self.assertEqual(self.lns.weight, 1)
        self.assertEqual(self.lns.zeta, 0.1)
        self.assertEqual(self.lns.sigma, 2)

    def test_pickling(self):
        p = pickle.dumps(self.lns, -1)
        up = pickle.loads(p)
        self.assertEqual(up.beg, 0)
        self.assertEqual(up.end, 1)
        self.assertEqual(up.weight, 1)
        self.assertEqual(up.zeta, 0.1)
        self.assertEqual(up.sigma, 2)

    def test_call(self):
        m = call_Sregion.call(self.lns, 191)
        self.assertTrue(m.s >= 0.0)
        lns = fwdpy11.LogNormalS(0, 1, 1, 0.1, 2, scaling=-1)
        m = call_Sregion.call(lns, 191)
        self.assertTrue(m.s <= 0.0)


class TestMultivariateLogNormalS(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.lns = fwdpy11.LogNormalS.mv(0, 1, 1, scaling=3.0, coupled=False)

    def test_create(self):
        self.assertEqual(self.lns.beg, 0)
        self.assertEqual(self.lns.end, 1)
        self.assertEqual(self.lns.weight, 1)
        self.assertEqual(self.lns.coupled, False)
        self.assertTrue(self.lns.zeta is None)
        self.assertTrue(self.lns.sigma is None)

    def test_pickling(self):
        p = pickle.dumps(self.lns, -1)
        up = pickle.loads(p)
        self.assertEqual(up.beg, 0)
        self.assertEqual(up.end, 1)
        self.assertEqual(up.weight, 1)
        self.assertEqual(up.coupled, False)
        self.assertEqual(up.scaling, 3.0)
        self.assertTrue(up.zeta is None)
        self.assertTrue(up.sigma is None)

    def test_call(self):
        mv = fwdpy11.mvDES(self.lns, np.zeros(2), np.identity(2))
        m = call_Sregion.call(mv, 191)
        self.assertTrue(np.all(m.esizes >= 0.0))
        lns = fwdpy11.LogNormalS.mv(0, 1, 1, scaling=-1)
        mv = fwdpy11.mvDES(lns, np.zeros(2), np.identity(2))
        m = call_Sregion.call(mv, 191)
        self.assertTrue(np.all(m.esizes <= 0.0))


class testPoissonInterval(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.pi = fwdpy11.PoissonInterval(0, 1, 1e-3)

    def test_pickling(self):
        p = pickle.dumps(self.pi)
        up = pickle.loads(p)
        self.assertEqual(up.beg, self.pi.beg)
        self.assertEqual(up.end, self.pi.end)
        self.assertEqual(up.mean, self.pi.mean)


class testBinomialInterval(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.bi = fwdpy11.BinomialInterval(0, 1, 1e-3)

    def test_pickling(self):
        p = pickle.dumps(self.bi)
        up = pickle.loads(p)
        self.assertEqual(up.beg, self.bi.beg)
        self.assertEqual(up.end, self.bi.end)
        self.assertEqual(up.probability, self.bi.probability)


class testFixedCrossovers(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.pi = fwdpy11.FixedCrossovers(0, 1, 10)

    def test_pickling(self):
        p = pickle.dumps(self.pi)
        up = pickle.loads(p)
        self.assertEqual(up.beg, self.pi.beg)
        self.assertEqual(up.end, self.pi.end)
        self.assertEqual(up.num_xovers, self.pi.num_xovers)


class testPoissonPoint(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.pi = fwdpy11.PoissonPoint(1.0, 0.5)

    def test_pickling(self):
        p = pickle.dumps(self.pi)
        up = pickle.loads(p)
        self.assertEqual(up.position, self.pi.position)
        self.assertEqual(up.mean, self.pi.mean)


class testBinomialPoint(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.pi = fwdpy11.BinomialPoint(1.0, 0.5)

    def test_pickling(self):
        p = pickle.dumps(self.pi)
        up = pickle.loads(p)
        self.assertEqual(up.position, self.pi.position)
        self.assertEqual(up.probability, self.pi.probability)


if __name__ == "__main__":
    unittest.main()
