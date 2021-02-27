import lzma
import pickle
import unittest

import fwdpy11


class TestOldBinaryFileFormats(unittest.TestCase):
    def test_reading_045(self):
        pop = fwdpy11.DiploidPopulation.load_from_file("tests/v045.bin")
        assert isinstance(pop, fwdpy11.DiploidPopulation) is True
        for i in pop.tables.mutations:
            self.assertEqual(
                pop.mutations[i.key].pos, pop.tables.sites[i.site].position
            )

if __name__ == "__main__":
    unittest.main()
