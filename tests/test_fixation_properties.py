# Test behavior of fixations in simulations.
# As of 0.1.3a1, fixations are to be sorted by position
# for all types of simulations, allowing for easier lookup.

import unittest
import os
from quick_pops import quick_nonneutral_slocus
from quick_pops import quick_mlocus_qtrait_change_optimum


class testFixationsAreSorted(unittest.TestCase):
    def test_fixations_SlocusPop(self):
        from fwdpy11.regions import ExpS
        pop = quick_nonneutral_slocus(N=1000,
                                      simlen=10000,
                                      dfe=ExpS(0, 1, 1, 0.05))
        self.assertTrue(len(pop.fixations) > 0)
        self.assertEqual(len(pop.fixations), len(pop.fixation_times))
        fpos = [i.pos for i in pop.fixations]
        self.assertTrue(sorted(fpos))
        non_neutral_fixations = [i.key for i in pop.fixations if i.neutral is False]
        # If this test fails, we have sim parameters
        # that are not useful for testing:
        self.assertTrue(len(non_neutral_fixations) > 0)
        for ni in non_neutral_fixations:
            self.assertFalse(any(i.key == ni for i in pop.mutations))


class testFixationsAreSortedQtraitSim(unittest.TestCase):
    """
    In sims of quantitative traits, selected mutations are
    retained, as they still affect mean trait value, and thus
    mean distance from optimum, and thus mean fitness.
    However, they are entered into pop.fixations so that
    their fixation times may be recorded.
    """

    @unittest.skipIf("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
                     "Skipping this test on Travis CI.")
    def test_sorted_fixations_MlocusPopQtrait(self):
        pop = quick_mlocus_qtrait_change_optimum(N=1000, simlen=10000)
        self.assertTrue(len(pop.fixations) > 0)
        self.assertEqual(len(pop.fixations), len(pop.fixation_times))
        fpos = [i.pos for i in pop.fixations]
        self.assertTrue(sorted(fpos))
        non_neutral_fixations = [i.key for i in pop.fixations if i.neutral is False]
        # If this test fails, we have sim parameters
        # that are not useful for testing:
        self.assertTrue(len(non_neutral_fixations) > 0)
        for ni in non_neutral_fixations:
            # These non-neutral fixations
            # must still be found in the pop.
            self.assertTrue(any(i.key == ni for i in pop.mutations))


if __name__ == "__main__":
    unittest.main()
