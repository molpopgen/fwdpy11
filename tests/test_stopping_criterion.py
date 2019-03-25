import unittest
import numpy as np
import fwdpy11


class test_stopping_criterion_DiploidPopulation(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.pop = fwdpy11.DiploidPopulation(1000, 1.0)
        p = {'nregions': [],  # No neutral mutations -- add them later!
             'gvalue': fwdpy11.Additive(2.0),
             'sregions': [fwdpy11.ExpS(0, 1, 1, -0.1)],
             'recregions': [fwdpy11.Region(0, 1, 1)],
             'rates': (0.0, 1e-3, 1e-3),
             # Keep mutations at frequency 1 in the pop if they affect fitness.
             'prune_selected': False,
             'demography':  np.array([1000]*10000, dtype=np.uint32)
             }
        self.params = fwdpy11.ModelParams(**p)

    def test_stop(self):
        rng = fwdpy11.GSLrng(42)

        def generation(pop, simplified):
            if pop.generation >= 50:
                return True
            return False

        fwdpy11.evolvets(
            rng, self.pop, self.params, 100,
            stopping_criterion=generation)
        self.assertEqual(self.pop.generation, 50)


if __name__ == "__main__":
    unittest.main()
