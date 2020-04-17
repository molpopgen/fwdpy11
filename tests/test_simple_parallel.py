# concurrent.futures is Python 3 only
import concurrent.futures as cf
import unittest

import numpy as np

import fwdpy11 as fp11
import fwdpy11.ezparams as fp11ez


def evolve_and_return(args):
    """
    This function runs our simulation.
    The input arguments come in a tuple,
    which is required by many of Python's
    functions for execution in separate processes.

    For this function, the arguments are the population
    size and a random number seed.
    """
    from fwdpy11 import Multiplicative

    N, seed = args
    # Construct as single-deme object
    # with N diploids
    pop = fp11.DiploidPopulation(N)
    theta = 100.0
    # Initialize a random number generator
    rng = fp11.GSLrng(seed)
    p = fp11ez.mslike(
        pop,
        simlen=100,
        rates=(theta / float(4 * pop.N), 1e-3, theta / float(4 * pop.N)),
    )
    p["gvalue"] = Multiplicative(2.0)
    params = fp11.ModelParams(**p)
    fp11.evolve_genomes(rng, pop, params)
    # The population is picklable, and so
    # we can return it from another process
    return pop


class testSimpleParallel(unittest.TestCase):
    """
    This would be your __name__=="__main__"
    in a script
    """

    def test_simple_parallel(self):
        np.random.seed(101)
        # Generate a list of arguments for our processes.
        # We generation 10 random number seeds
        args = [(1000, seed) for seed in np.random.randint(0, 42000000, 10)]
        # Use a process pool with a max of 10 workers
        with cf.ProcessPoolExecutor(10) as pool:
            # Run our simulations and get the
            # result back, which will be
            # the population
            for res in pool.map(evolve_and_return, args):
                self.assertEqual(res.N, 1000)
                self.assertEqual(res.generation, 100)


if __name__ == "__main__":
    unittest.main()
