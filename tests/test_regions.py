# Unit tests for fwdpy11.regions

import fwdpy11
import numpy as np
import unittest


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


class testConstantS(unittest.TestCase):
    def testBadValue(self):
        with self.assertRaises(ValueError):
            fwdpy11.ConstantS(0, 1, 1, np.nan)

    def testBadDominance(self):
        with self.assertRaises(ValueError):
            fwdpy11.ConstantS(0, 1, 1, 1, np.nan)


class testMutationRegions(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.nregions = [fwdpy11.Region(0, 1, 1)]
        self.sregions = [fwdpy11.ConstantS(
            0, 1, 1, -0.1, 0.25), fwdpy11.ExpS(0, 1, 0.001, 0.01, 1.0)]

    def test_create(self):
        mr = fwdpy11.MutationRegions.create(0.5, self.nregions, self.sregions)  # NOQA


class testMlocusMutationRegions(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.nregions = [fwdpy11.Region(0, 1, 1)]
        self.sregions = [fwdpy11.ConstantS(
            0, 1, 1, -0.1, 0.25), fwdpy11.ExpS(0, 1, 0.001, 0.01, 1.0)]

    def test_create(self):
        mr = [fwdpy11.MutationRegions.create(
            0.5, self.nregions, self.sregions)] * 5
        mmr = fwdpy11.MlocusMutationRegions()
        for i in mr:
            mmr.append(i)
        self.assertEqual(len(mmr), 5)


class testRecombinationRegions(unittest.TestCase):
    def test_create(self):
        rr = fwdpy11.RecombinationRegions(1e-3,
                                          [fwdpy11.Region(0, 1, 1),
                                              fwdpy11.Region(1, 2, 1)])
        self.assertEqual(len(rr.weights), 2)

    def test_bad_input(self):
        with self.assertRaises(TypeError):
            rr = fwdpy11.RecombinationRegions(1e-3,
                                              [fwdpy11.ExpS(0, 1, 1, -0.2), fwdpy11.GaussianS(1, 2, 1, 0.25)])


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
