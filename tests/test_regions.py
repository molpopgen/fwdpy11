# Unit tests for fwdpy11.regions

import fwdpy11
import unittest


class test_Region_repr(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.nregion = fwdpy11.Region(0, 1, 1)
        self.constants = fwdpy11.ConstantS(0, 1, 1, -0.1, 2.0)
        self.uniforms = fwdpy11.UniformS(0, 1, 1, 1, 2, 0)
        self.exps = fwdpy11.ExpS(0, 1, 1, 0.25)
        self.gaussians = fwdpy11.GaussianS(0, 1, 1, 0.25)
        self.gammas = fwdpy11.GammaS(0, 1, 1, 0.124, 2.0)

    def test_nregion(self):
        x = eval('fwdpy11.' + repr(self.nregion))
    def test_constants(self):
        x = eval('fwdpy11.' + repr(self.constants))
    def test_uniforms(self):
        x = eval('fwdpy11.' + repr(self.uniforms))
    def test_exps(self):
        x = eval('fwdpy11.' + repr(self.exps))
    def test_gaussians(self):
        x = eval('fwdpy11.' + repr(self.gaussians))
    def test_gammas(self):
        x = eval('fwdpy11.' + repr(self.gammas))

if __name__ == "__main__":
    unittest.main()
