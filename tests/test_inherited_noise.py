import unittest

import fwdpy11
import fwdpy11.custom_genetic_value_decorators
import numpy as np


@fwdpy11.custom_genetic_value_decorators.genetic_value_noise_default_clone
class InheritedNoise(fwdpy11.GeneticValueNoise):
    def __init__(self):
        self.generation = 0
        fwdpy11.GeneticValueNoise.__init__(self)

    def __call__(self, data) -> float:
        return data.parent1_metadata.e + self.generation

    def update(self, pop):
        self.generation = pop.generation


class TestInheritNoise(unittest.TestCase):
    """
    The class InheritedNoise passes the
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
        self.pop = fwdpy11.DiploidPopulation([100, 100, 100], 1.0)
        self.pdict = {
            "nregions": [],
            "sregions": [],
            "recregions": [],
            "rates": (0, 0, 0),
            "demography": fwdpy11.ForwardDemesGraph.tubes(
                [100] * 3, 1),
            "simlen": 3,
            "gvalue": fwdpy11.Additive(
                2.0, fwdpy11.GSS(optimum=0.0, VS=1.0), InheritedNoise()
            ),
        }

    def test_noise_values_tree_sequences(self):
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)
        self.assertEqual(self.pop.generation, 3)
        self.assertTrue(
            all([i.e == 6.0 for i in self.pop.diploid_metadata]) is True)


if __name__ == "__main__":
    unittest.main()
