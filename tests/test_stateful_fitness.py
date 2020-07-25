import math
import pickle
import unittest

import numpy as np

import fwdpy11 as fp11
import fwdpy11.ezparams
import snowdrift


class SamplePhenotypes(object):
    """
    Temporal sampler checks that one can hook
    a stateful fitness model to a sampler
    and access its data and that the data
    are as expected.
    """

    def __init__(self, f, slope, p0):
        self.f = f
        self.slope = slope
        self.p0 = p0
        self.sig0 = (1.0 / slope) * math.log(self.p0 / (1.0 - self.p0))

    def __call__(self, pop, sampler):
        for i in range(pop.N):
            w = 0.0
            for m in pop.haploid_genomes[pop.diploids[i].first].smutations:
                w += pop.mutations[m].s
            for m in pop.haploid_genomes[pop.diploids[i].second].smutations:
                w += pop.mutations[m].s
            w = 1.0 / (1.0 + math.exp(-self.slope * (w + self.sig0)))
            assert np.isclose(w, self.f.phenotypes[i])


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
    pop = fp11.DiploidPopulation(N, 1.0)
    # Initialize a random number generator
    rng = fp11.GSLrng(seed)
    p = {
        "sregions": [fp11.ExpS(0, 1, 1, -0.1, 1.0)],
        "recregions": [fp11.Region(0, 1, 1)],
        "nregions": [],
        "gvalue": snowdrift.DiploidSnowdrift(0.2, -0.2, 1, -2, 1, 0.3),
        # evolve for 100 generations so that unit tests are
        # fast
        "demography": fwdpy11.DiscreteDemography(),
        "simlen": 20,
        "rates": (0.0, 0.0025, 0.001),
        "prune_selected": False,
    }
    params = fwdpy11.ModelParams(**p)
    sampler = SamplePhenotypes(params.gvalue, 1, 0.3)
    fp11.evolvets(rng, pop, params, 100, sampler)
    # return our pop
    return pop


class testSnowdrift(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.f = snowdrift.DiploidSnowdrift(1, -1, 0.1, 0.2, 1, 0.3)

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
        self.assertEqual(up.phenotypes, self.f.phenotypes)

    def test_evolve(self):
        p = evolve_snowdrift((1000, 42))
        p


if __name__ == "__main__":
    unittest.main()
