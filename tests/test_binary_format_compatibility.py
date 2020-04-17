import lzma
import pickle
import unittest

import fwdpy11


class TestOldBinaryFileFormats(unittest.TestCase):
    def test_reading_045(self):
        pop = fwdpy11.DiploidPopulation.load_from_file("tests/v045.bin")
        for i in pop.tables.mutations:
            self.assertEqual(
                pop.mutations[i.key].pos, pop.tables.sites[i.site].position
            )

    def test_unpickling_045(self):
        with lzma.open("tests/v045.lzma", "rb") as f:
            pop = pickle.load(f)
        for i in pop.tables.mutations:
            self.assertEqual(
                pop.mutations[i.key].pos, pop.tables.sites[i.site].position
            )


if __name__ == "__main__":
    unittest.main()
