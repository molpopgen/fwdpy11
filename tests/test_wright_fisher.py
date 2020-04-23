import unittest

import numpy as np

import fwdpy11


class GenerationRecorder(object):
    def __init__(self):
        self.generations = []

    def __call__(self, pop):
        self.generations.append(pop.generation)


class testWFevolve(unittest.TestCase):
    @classmethod
    def setUpClass(self):

        self.pop = fwdpy11.DiploidPopulation(1000)
        self.rng = fwdpy11.GSLrng(42)
        self.recorder = GenerationRecorder()
        self.p = {
            "rates": (1e-3, 1e-3, 1e-3),
            "popsizes": np.array([1000] * 100, dtype=np.uint32),
            "simlen": 100,
            "nregions": [fwdpy11.Region(0, 1, 1)],
            "sregions": [fwdpy11.ExpS(0, 1, 1, -1e-2)],
            "recregions": [fwdpy11.Region(0, 1, 1)],
            "gvalue": fwdpy11.Multiplicative(2.0),
        }

    def testEvolve(self):

        params = fwdpy11.ModelParams(**self.p)
        fwdpy11.evolve_genomes(self.rng, self.pop, params, self.recorder)
        self.assertEqual(self.recorder.generations, [i + 1 for i in range(100)])
        self.p["popsizes"] = np.array([self.pop.N] * 24, dtype=np.uint32)
        self.p["simlen"] = len(self.p["popsizes"])
        params = fwdpy11.ModelParams(**self.p)
        fwdpy11.evolve_genomes(self.rng, self.pop, params, self.recorder)
        self.assertEqual(self.recorder.generations, [i + 1 for i in range(124)])


if __name__ == "__main__":
    unittest.main()
