import unittest

import fwdpy11
import numpy as np


class test_stopping_criterion_DiploidPopulation(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.pop = fwdpy11.DiploidPopulation(1000, 1.0)
        p = {
            "nregions": [],  # No neutral mutations -- add them later!
            "gvalue": fwdpy11.Additive(2.0),
            "sregions": [fwdpy11.ExpS(0, 1, 1, -0.1)],
            "recregions": [fwdpy11.PoissonInterval(0, 1, 1e-3)],
            "rates": (0.0, 1e-3, None),
            # Keep mutations at frequency 1 in the pop if they affect fitness.
            "prune_selected": False,
            "demography": None,
            "simlen": 10 * self.pop.N,
        }
        self.params = fwdpy11.ModelParams(**p)

    def test_stop(self):
        rng = fwdpy11.GSLrng(42)

        def generation(pop, simplified):
            if pop.generation >= 50:
                return True
            return False

        fwdpy11.evolvets(rng, self.pop, self.params, 100, stopping_criterion=generation)
        self.assertEqual(self.pop.generation, 50)


def test_stopping_criterion_with_ancient_samples():
    """
    This triggers github issue 844
    """

    def preserve(pop, sampler):
        if pop.generation == 10:
            sampler.assign(np.arange(pop.N, dtype=np.uint32))

    def stopper(pop, _):
        if pop.generation == 10:
            return True
        return False

    pdict = {
        "sregions": [fwdpy11.ExpS(0, 1, 1, -0.05, 1)],
        "gvalue": fwdpy11.Multiplicative(2.0),
        "rates": (0, 0.1, 0),
        "simlen": 20,
    }
    params = fwdpy11.ModelParams(**pdict)
    pop = fwdpy11.DiploidPopulation(100, 1.0)
    rng = fwdpy11.GSLrng(90210)

    fwdpy11.evolvets(
        rng, pop, params, 100, recorder=preserve, stopping_criterion=stopper
    )


if __name__ == "__main__":
    unittest.main()
