# Test behavior of fixations in simulations.
# As of 0.1.3a1, fixations are to be sorted by position
# for all types of simulations, allowing for easier lookup.

import unittest
import os
import fwdpy11
import fwdpy11.util
import fwdpy11.wright_fisher_qtrait
import numpy as np
from quick_pops import quick_nonneutral_slocus
from quick_pops import quick_slocus_qtrait_pop_params
from quick_pops import quick_mlocus_qtrait_change_optimum


class testFixationsAreSorted(unittest.TestCase):
    def test_fixations_SlocusPop(self):
        from fwdpy11 import ExpS
        pop = quick_nonneutral_slocus(N=1000,
                                      simlen=10000,
                                      dfe=ExpS(0, 1, 1, 0.05))
        self.assertTrue(len(pop.fixations) > 0)
        self.assertEqual(len(pop.fixations), len(pop.fixation_times))
        fpos = [i.pos for i in pop.fixations]
        self.assertTrue(sorted(fpos))
        non_neutral_fixations = [
            i.key for i in pop.fixations if i.neutral is False]
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
    def test_sorted_fixations_MlocusPopQtrait_with_pruning(self):
        pop = quick_mlocus_qtrait_change_optimum(N=1000, simlen=10000)
        self.assertTrue(len(pop.fixations) > 0)
        self.assertEqual(len(pop.fixations), len(pop.fixation_times))
        fpos = [i.pos for i in pop.fixations]
        self.assertTrue(sorted(fpos))
        non_neutral_fixations = [
            i.key for i in pop.fixations if i.neutral is False]
        # If this test fails, we have sim parameters
        # that are not useful for testing:
        self.assertTrue(len(non_neutral_fixations) > 0)
        for ni in non_neutral_fixations:
            # These non-neutral fixations
            # must still be found in the pop.
            self.assertTrue(any(i.key == ni for i in pop.mutations))
        neutral_fixations = [i.key for i in pop.fixations if i.neutral is True]
        self.assertTrue(len(neutral_fixations) > 0)
        for ni in neutral_fixations:
            self.assertTrue(any(i.key == ni for i in pop.mutations) is False)

    @unittest.skipIf("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
                     "Skipping this test on Travis CI.")
    def test_sorted_fixations_MlocusPopQtrait_without_pruning(self):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pop = quick_mlocus_qtrait_change_optimum(N=1000, simlen=10000, prune_selected = True)
            self.assertTrue(len(pop.fixations) > 0)
            self.assertEqual(len(pop.fixations), len(pop.fixation_times))
            fpos = [i.pos for i in pop.fixations]
            self.assertTrue(sorted(fpos))
            non_neutral_fixations = [
                i.key for i in pop.fixations if i.neutral is False]
            # If this test fails, we have sim parameters
            # that are not useful for testing:
            self.assertTrue(len(non_neutral_fixations) > 0)
            for ni in non_neutral_fixations:
                # These non-neutral fixations
                # must not be found in the pop.
                self.assertTrue(any(i.key == ni for i in pop.mutations) is False)
            neutral_fixations = [i.key for i in pop.fixations if i.neutral is True]
            self.assertTrue(len(neutral_fixations) > 0)
            for ni in neutral_fixations:
                self.assertTrue(any(i.key == ni for i in pop.mutations) is False)

class tests_MultipleFixationsAtPosition(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        from fwdpy11.model_params import SlocusParamsQ
        self.pop, self.pdict = quick_slocus_qtrait_pop_params()
        self.params = SlocusParamsQ(**self.pdict)
        self.params.demography = np.array([self.pop.N], dtype=np.uint32)
        self.params.mutrate_n = 0.0  # disallow neutral mutations
        self.rng = fwdpy11.GSLrng(42)
        fwdpy11.util.add_mutation(
            self.rng, self.pop, 2 * self.pop.N, (0.1, 0.0, 0.0), 0)

    def test_two_neutral_fixations(self):
        # Evolve pop for 1 gen.
        # This will push mutation into fixations:
        fwdpy11.wright_fisher_qtrait.evolve(self.rng, self.pop, self.params)
        self.assertEqual(len(self.pop.fixations), 1)
        # Add back the exact same mutation:
        fwdpy11.util.add_mutation(
            self.rng, self.pop, 2 * self.pop.N, (0.1, 0.0, 0.0), 0)
        # Evolve again for 1 gen:
        fwdpy11.wright_fisher_qtrait.evolve(self.rng, self.pop, self.params)
        self.assertEqual(len(self.pop.fixations), 2)
        self.assertEqual(list(self.pop.fixation_times), [1, 2])


if __name__ == "__main__":
    unittest.main()
