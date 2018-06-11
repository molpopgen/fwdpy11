import unittest
import fwdpy11.genetic_values
import numpy as np


class TestSlocusAdditive(unittest.TestCase):
    """
    API and behavior tests
    """

    def testConstructionWithNumpyFloat(self):
        x = fwdpy11.genetic_values.SlocusAdditive(np.float(1.0))
        self.assertEqual(x.scaling, 1.0)

    def testConstructWithNaN(self):
        with self.assertRaises(ValueError):
            x = fwdpy11.genetic_values.SlocusAdditive(np.nan)


if __name__ == "__main__":
    unittest.main()
