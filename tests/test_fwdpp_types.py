import unittest
from quick_pops import quick_nonneutral_slocus


class testWriteMutationLabel(unittest.TestCase):
    def test_write_label(self):
        pop = quick_nonneutral_slocus()
        for i in pop.mutations:
            if i.neutral is False:
                i.label = 3
            else:
                i.label = 17
        for i in pop.mutations:
            if i.neutral is False:
                self.assertEqual(i.label, 3)
            else:
                self.assertEqual(i.label, 17)


if __name__ == "__main__":
    unittest.main()
