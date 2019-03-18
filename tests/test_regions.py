# Unit tests for fwdpy11.regions

import fwdpy11
import numpy as np
import unittest
import pickle

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
        self.assertEqual(self.r.b, up.b)
        self.assertEqual(self.r.e, up.e)
        self.assertEqual(self.r.w, up.w)
        self.assertEqual(self.r.c, up.c)
        self.assertEqual(self.r.l, up.l)


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
        self.assertEqual(self.r.b, up.b)
        self.assertEqual(self.r.e, up.e)
        self.assertEqual(self.r.w, up.w)
        self.assertEqual(self.r.c, up.c)
        self.assertEqual(self.r.l, up.l)
        self.assertEqual(up.mean, -0.2)
        self.assertEqual(up.h, DOM)


class testGammaS(unittest.TestCase):
    def testBadMean(self):
        with self.assertRaises(ValueError):
            fwdpy11.GammaS(0, 1, 1, mean=np.nan, shape=1.0)

    def testBadShape(self):
        with self.assertRaises(ValueError):
            fwdpy11.GammaS(0, 1, 1, mean=1., shape=np.nan)

    def testBadDominance(self):
        with self.assertRaises(ValueError):
            fwdpy11.GammaS(0, 1, 1, mean=1., shape=1., h=np.nan)


class test_PickleGammaS(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.r = fwdpy11.GammaS(BEG, END, WEIGHT,
                                -0.2, 0.5, DOM, COUPLED, LABEL)

    def test_pickling(self):
        p = pickle.dumps(self.r, -1)
        up = pickle.loads(p)
        self.assertEqual(self.r.b, up.b)
        self.assertEqual(self.r.e, up.e)
        self.assertEqual(self.r.w, up.w)
        self.assertEqual(self.r.c, up.c)
        self.assertEqual(self.r.l, up.l)
        self.assertEqual(up.mean, -0.2)
        self.assertEqual(up.shape, 0.5)
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
            fwdpy11.UniformS(0, 1, 1, lo=1., hi=1., h=np.nan)


class test_PickleUniformS(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.r = fwdpy11.UniformS(BEG, END, WEIGHT,
                                  -0.2, 0.5, DOM, COUPLED, LABEL)

    def test_pickling(self):
        p = pickle.dumps(self.r, -1)
        up = pickle.loads(p)
        self.assertEqual(self.r.b, up.b)
        self.assertEqual(self.r.e, up.e)
        self.assertEqual(self.r.w, up.w)
        self.assertEqual(self.r.c, up.c)
        self.assertEqual(self.r.l, up.l)
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
        self.r = fwdpy11.GaussianS(BEG, END, WEIGHT,
                                   0.5, DOM, COUPLED, LABEL)

    def test_pickling(self):
        p = pickle.dumps(self.r, -1)
        up = pickle.loads(p)
        self.assertEqual(self.r.b, up.b)
        self.assertEqual(self.r.e, up.e)
        self.assertEqual(self.r.w, up.w)
        self.assertEqual(self.r.c, up.c)
        self.assertEqual(self.r.l, up.l)
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
        self.r = fwdpy11.ConstantS(BEG, END, WEIGHT,
                                   0.5, DOM, COUPLED, LABEL)

    def test_pickling(self):
        p = pickle.dumps(self.r, -1)
        up = pickle.loads(p)
        self.assertEqual(self.r.b, up.b)
        self.assertEqual(self.r.e, up.e)
        self.assertEqual(self.r.w, up.w)
        self.assertEqual(self.r.c, up.c)
        self.assertEqual(self.r.l, up.l)
        self.assertEqual(up.s, 0.5)
        self.assertEqual(up.h, DOM)


class testMutationRegions(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.nregions = [fwdpy11.Region(0, 1, 1)]
        self.sregions = [fwdpy11.ConstantS(
            0, 1, 1, -0.1, 0.25), fwdpy11.ExpS(0, 1, 0.001, 0.01, 1.0)]

    def test_create(self):
        mr = fwdpy11.MutationRegions.create(0.5, self.nregions, self.sregions)  # NOQA


class testMultivariateGaussianEffects(unittest.TestCase):
    def testConstruction(self):
        try:
            fwdpy11.MultivariateGaussianEffects(
                0, 1, 1, np.identity(2))
        except Exception:
            self.fail("Unexpected exception during object construction")

    def testInvalidMatrixShape(self):
        # Input matrix has wrong shape:
        m = np.array([-1.0] * 4)
        with self.assertRaises(ValueError):
            fwdpy11.MultivariateGaussianEffects(
                0, 1, 1, m)

    def testInvalidMatrixShape2(self):
        m = np.array([-1.0] * 6).reshape((3, 2))
        with self.assertRaises(ValueError):
            fwdpy11.MultivariateGaussianEffects(
                0, 1, 1, m)

    def testInvalidMatrix(self):
        # This is not a positive-definite covariance
        # matrix:
        m = np.array([-1.0] * 4).reshape((2, 2))
        with self.assertRaises(ValueError):
            fwdpy11.MultivariateGaussianEffects(
                0, 1, 1, m)

    def testMatrixWithNAN(self):
        m = np.identity(4).reshape((4, 4))
        m[0, 1] = np.nan
        with self.assertRaises(ValueError):
            fwdpy11.MultivariateGaussianEffects(0, 1, 1, m)


class test_PickleMultivariateGaussianEffects(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.r = fwdpy11.MultivariateGaussianEffects(
            BEG, END, WEIGHT, np.identity(2), -0.3, DOM, COUPLED, LABEL)

    def test_pickling(self):
        p = pickle.dumps(self.r, -1)
        up = pickle.loads(p)
        self.assertEqual(self.r, up)


class testRecombinationRegions(unittest.TestCase):
    def test_create(self):
        rr = fwdpy11.RecombinationRegions(1e-3,
                                          [fwdpy11.Region(0, 1, 1),
                                              fwdpy11.Region(1, 2, 1)])
        self.assertEqual(len(rr.weights), 2)

    def test_bad_input(self):
        with self.assertRaises(TypeError):
            fwdpy11.RecombinationRegions(1e-3,
                                         [fwdpy11.ExpS(0, 1, 1, -0.2),
                                          fwdpy11.GaussianS(1, 2, 1, 0.25)])


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


class testFixedCrossovers(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.pi = fwdpy11.FixedCrossovers(0, 1, 10)

    def test_pickling(self):
        p = pickle.dumps(self.pi)
        up = pickle.loads(p)
        self.assertEqual(up.beg, self.pi.beg)
        self.assertEqual(up.end, self.pi.end)
        self.assertEqual(up.nxovers, self.pi.nxovers)


class testPoissonPoint(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.pi = fwdpy11.PoissonPoint(1., 0.5)

    def test_pickling(self):
        p = pickle.dumps(self.pi)
        up = pickle.loads(p)
        self.assertEqual(up.position, self.pi.position)
        self.assertEqual(up.mean, self.pi.mean)


class testBinomialPoint(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.pi = fwdpy11.BinomialPoint(1., 0.5)

    def test_pickling(self):
        p = pickle.dumps(self.pi)
        up = pickle.loads(p)
        self.assertEqual(up.position, self.pi.position)
        self.assertEqual(up.probability, self.pi.probability)

# class test_Region_repr(unittest.TestCase):
#     @classmethod
#     def setUpClass(self):
#         self.nregion = fwdpy11.Region(0, 1, 1)
#         self.constants = fwdpy11.ConstantS(0, 1, 1, -0.1, 2.0)
#         self.uniforms = fwdpy11.UniformS(0, 1, 1, 1, 2, 0)
#         self.exps = fwdpy11.ExpS(0, 1, 1, 0.25)
#         self.gaussians = fwdpy11.GaussianS(0, 1, 1, 0.25)
#         self.gammas = fwdpy11.GammaS(0, 1, 1, 0.124, 2.0)
#
#     def test_nregion(self):
#         x = eval('fwdpy11.' + repr(self.nregion))
#
#     def test_constants(self):
#         x = eval('fwdpy11.' + repr(self.constants))
#
#     def test_uniforms(self):
#         x = eval('fwdpy11.' + repr(self.uniforms))
#
#     def test_exps(self):
#         x = eval('fwdpy11.' + repr(self.exps))
#
#     def test_gaussians(self):
#         x = eval('fwdpy11.' + repr(self.gaussians))
#
#     def test_gammas(self):
#         x = eval('fwdpy11.' + repr(self.gammas))
#
#
# class test_Scaling(unittest.TestCase):
#     def test_ConstantSwithoutScaling(self):
#         """
#         Here, the scaling is 1, and the DFE is constant.
#         Thus, all "s" values must be equal to 10.0
#         """
#         x = fwdpy11.ConstantS(0, 1, 1, s=10, scaling=1)
#         c = x.callback()
#         rng = fwdpy11.GSLrng(42)
#         sh = c(rng)
#         self.assertEqual(sh[0], 10.0)
#
#     def test_ConstantSwithScaling(self):
#         """
#         Here, the scaling is 100, and the DFE is constant.
#         Thus, all "s" values must be equal to 10.0/100.0,
#         or 0.1.
#         """
#         x = fwdpy11.ConstantS(0, 1, 1, s=10, scaling=100)
#         c = x.callback()
#         rng = fwdpy11.GSLrng(42)
#         sh = c(rng)
#         self.assertEqual(sh[0], 0.1)


if __name__ == "__main__":
    unittest.main()
