import unittest

import fwdpy11

from utils import make_path


class TestOldBinaryFileFormats(unittest.TestCase):
    def test_reading_045(self):
        path_to_file = make_path("v045.bin")
        pop = fwdpy11.DiploidPopulation.load_from_file(path_to_file)
        assert isinstance(pop, fwdpy11.DiploidPopulation) is True
        for i in pop.tables.mutations:
            self.assertEqual(
                pop.mutations[i.key].pos, pop.tables.sites[i.site].position
            )


if __name__ == "__main__":
    unittest.main()
