import unittest
import fwdpy11 as fp11
import fwdpy11.multilocus as fp11m
import numpy as np

class testRecombination(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.rng = fp11.GSLrng(42)

    def testMakeBinomial(self):
        v = [0.5] * 5
        x = fp11m.binomial_rec(v)

    def testMakePoisson(self):
        v = [1e-3] * 5
        x = fp11m.poisson_rec(v)


if __name__ == "__main__":
    unittest.main()

