import unittest
import fwdpy11
from quick_pops import quick_nonneutral_slocus


class test_mutation_labels(unittest.TestCase):
    def test_label_selected_mutations(self):
        pop = quick_nonneutral_slocus(
            dfe=fwdpy11.ExpS(0, 1, 1, -0.1, 1.0, label=2))
        selected = 0
        for i in pop.mutations:
            if i.neutral is True:
                self.assertEqual(i.label, 0)
            else:
                self.assertEqual(i.label, 2)
                selected += 1
        self.assertTrue(selected > 0) #Else test is useless


if __name__ == "__main__":
    unittest.main()
