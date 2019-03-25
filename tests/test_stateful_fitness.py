import unittest
import pickle
import numpy as np
import fwdpy11 as fp11
import fwdpy11
import fwdpy11.ezparams
import snowdrift


class SamplePhenotypes(object):
    """
    Temporal sampler checks that one can hook
    a stateful fitness model to a sampler
    and access its data and that the data
    are as expected.
    """

    def __init__(self, f):
        self.f = f
        # self.a = fwdpy11.genetic_values.DiploidAdditive(2.0)

    def __call__(self, pop):
        pass
        # for i in range(pop.N):
        #     w = self.a(pop.diploids[i], pop)
        #     assert(w == self.f.phenotypes[i])


def evolve_snowdrift(args):
    """
    We write the function taking a tuple
    out of habit, simplifying later
    integration with multiprocessing or
    concurrent.futures.
    """
    N, seed = args
    # Construct as single-deme object
    # with N diploids
    pop = fp11.DiploidPopulation(N)
    # Initialize a random number generator
    rng = fp11.GSLrng(seed)
    p = {'sregions': [fp11.ExpS(0, 1, 1, -0.1, 1.0)],
         'recregions': [fp11.Region(0, 1, 1)],
         'nregions': [],
         'gvalue': snowdrift.DiploidSnowdrift(0.2, -0.2, 1, -2),
         # evolve for 100 generations so that unit tests are
         # fast
         'demography': np.array([N] * 100, dtype=np.uint32),
         'rates': (0.0, 0.0025, 0.001),
         'prune_selected': False
         }
    params = fwdpy11.ModelParams(**p)
    sampler = SamplePhenotypes(params.gvalue)
    fp11.evolve_genomes(rng, pop, params, sampler)
    # return our pop
    return pop


class testSnowdrift(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.f = snowdrift.DiploidSnowdrift(1, -1, 0.1, 0.2)

    def testShape(self):
        s = self.f.shape
        self.assertEqual(len(s), 1)
        self.assertEqual(s[0], 1)

    def test_genetic_values(self):
        self.assertEqual(len(self.f.genetic_values), 1)

    def testPickle(self):
        self.f.phenotypes = [1, 2, 3, 4]
        p = pickle.dumps(self.f, -1)
        up = pickle.loads(p)
        self.assertEqual(up.phenotypes, self.f. phenotypes)

    def test_evolve(self):
        p = evolve_snowdrift((1000, 42))
        p


if __name__ == "__main__":
    unittest.main()
