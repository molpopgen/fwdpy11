import unittest
from fwdpy11 import GSLrng
from fwdpy11.wright_fisher_qtrait import evolve, GaussianNoise, GSS
from fwdpy11.trait_values import SlocusAdditiveTrait
from fwdpy11.multilocus import MultiLocusGeneticValue
from quick_pops import quick_mlocus_qtrait_pop_params


class GaussianNoiseUpdater(GaussianNoise):
    """
    This noise updater inherits from the
    standard fwdpy11 object.  It defines an
    updater function to record how many
    times it has been called.
    """

    def __init__(self, rng, sd, mean=0.0):
        self.ncalls = 0
        super(GaussianNoiseUpdater, self).__init__(rng, sd, mean)

    def update(self, generation):
        self.ncalls += 1

class GSSupdater(GSS):
    def __init__(self, VS, O):
        super(GSSupdater, self).__init__(VS,O)
        self.ncalls = 0

    def update(self, pop):
        self.ncalls += 1


class testNoiseUpdater(unittest.TestCase):
    def test_noise_updater(self):
        rng = GSLrng(42)
        pop, params = quick_mlocus_qtrait_pop_params()
        n = GaussianNoiseUpdater(rng, 0.1)
        trait2w = GSSupdater(1.0,0)
        params.noise = n
        params.trait2w = trait2w
        evolve(rng, pop, params)
        self.assertEqual(n.ncalls, pop.generation)
        self.assertEqual(trait2w.ncalls, pop.generation)


if __name__ == "__main__":
    unittest.main()
