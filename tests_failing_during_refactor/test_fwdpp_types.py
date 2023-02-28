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


class testGetMutationKey(unittest.TestCase):
    def test_get_popgenmut_key(self):
        pop = quick_nonneutral_slocus()
        for i, j in zip(pop.mutations, pop.mcounts):
            if j > 0:
                self.assertEqual(i.key, (i.pos, i.s, i.g))


if __name__ == "__main__":
    unittest.main()
