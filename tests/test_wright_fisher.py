import unittest
import fwdpy11 as fp11
import fwdpy11.wright_fisher as wf

class testEmptyPopSizes(unittest.TestCase):
    def testRaisesException(self):
        p=fp11.Spop(100)
        rng=fp11.GSLrng(42)
        with self.assertRaises(RuntimeError):
            wf.evolve(rng,p,[])
