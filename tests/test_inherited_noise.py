import unittest

import numpy as np

import fwdpy11
import inherit_noise


class TestInheritNoise(unittest.TestCase):
    """
    The class inherit_noise.TestInheritNoise passes the
    sum of all generations into the "e" field of
    offspring metadata.  Thus, if offspring metadata
    can see parental metadata, then we can test that the
    final alive individuals (in a sim of non-overlapping
    generations) all have e fields that are the sum of
    the generations that have occurred.
    """

    @classmethod
    def setUp(self):
        self.rng = fwdpy11.GSLrng(42)
        self.pop = fwdpy11.DiploidPopulation(100, 1.0)
        self.pdict = {
            "nregions": [],
            "sregions": [],
            "recregions": [],
            "rates": (0, 0, 0),
            "demography": fwdpy11.DiscreteDemography(
                set_deme_sizes=np.array([self.pop.N] * 3, dtype=np.uint32)
            ),
            "popsizes": np.array([self.pop.N] * 3, dtype=np.uint32),
            "simlen": 3,
            "gvalue": fwdpy11.Additive(
                2.0, fwdpy11.GSS(optimum=0.0, VS=1.0), inherit_noise.IneritedNoise()
            ),
        }

    def test_noise_values_tree_sequences(self):
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)
        self.assertEqual(self.pop.generation, 3)
        self.assertTrue(all([i.e == 6.0 for i in self.pop.diploid_metadata]) is True)

    def test_noise_values_without_tree_sequences(self):
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolve_genomes(self.rng, self.pop, params)
        self.assertEqual(self.pop.generation, 3)
        self.assertTrue(all([i.e == 6.0 for i in self.pop.diploid_metadata]) is True)


if __name__ == "__main__":
    unittest.main()
