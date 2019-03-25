import unittest
import fwdpy11 as fp11
import pyximport
import numpy as np
pyximport.install(setup_args={'include_dirs': np.get_include()})  # noqa
from MeanFitness import MeanFitness


class GenerationRecorder(object):
    def __init__(self):
        self.generations = []

    def __call__(self, pop):
        self.generations.append(pop.generation)


class testWFevolve(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        from fwdpy11 import ModelParams
        from fwdpy11 import Multiplicative
        self.pop = fp11.DiploidPopulation(1000)
        self.rng = fp11.GSLrng(42)
        self.recorder = GenerationRecorder()
        self.cython_recorder = MeanFitness()
        self.p = ModelParams()
        self.p.rates = (1e-3, 1e-3, 1e-3)
        self.p.demography = np.array([1000] * 100, dtype=np.uint32)
        self.p.nregions = [fp11.Region(0, 1, 1)]
        self.p.sregions = [fp11.ExpS(0, 1, 1, -1e-2)]
        self.p.recregions = self.p.nregions
        self.p.gvalue = Multiplicative(2.0)

    def testEvolve(self):
        from fwdpy11 import evolve_genomes as evolve
        evolve(self.rng, self.pop, self.p, self.recorder)
        self.assertEqual(self.recorder.generations,
                         [i + 1 for i in range(100)])
        self.p.demography = np.array([self.pop.N] * 24)
        evolve(self.rng, self.pop, self.p, self.recorder)
        self.assertEqual(self.recorder.generations,
                         [i + 1 for i in range(124)])


class testCythonRecorder(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        from fwdpy11 import ModelParams
        from fwdpy11 import Multiplicative
        self.pop = fp11.DiploidPopulation(1000)
        self.rng = fp11.GSLrng(42)
        self.cython_recorder = MeanFitness()
        self.p = ModelParams()
        self.p.rates = (1e-3, 1e-3, 1e-3)
        self.p.demography = np.array([1000] * 100, dtype=np.uint32)
        self.p.nregions = [fp11.Region(0, 1, 1)]
        self.p.sregions = [fp11.ExpS(0, 1, 1, -1e-2)]
        self.p.recregions = self.p.nregions
        self.p.gvalue = Multiplicative(2.0)

    def testEvolve(self):
        from fwdpy11 import evolve_genomes as evolve
        evolve(self.rng, self.pop, self.p, self.cython_recorder)


if __name__ == "__main__":
    unittest.main()
