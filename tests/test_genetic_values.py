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
        self.assertEqual(x.is_fitness, True)

    def testConstructWithNaN(self):
        with self.assertRaises(ValueError):
            x = fwdpy11.genetic_values.SlocusAdditive(np.nan)

    def testConstructTrait(self):
        x = fwdpy11.genetic_values.SlocusAdditive(1.0, fwdpy11.genetic_values.AdditivePolicy.atrait,
                                                  fwdpy11.genetic_values.GSS(1.0, 0.0))
        self.assertEqual(x.scaling, 1.0)
        self.assertEqual(x.is_fitness, False)
        self.assertTrue(isinstance(x.gvalue_to_fitness,
                                   fwdpy11.genetic_values.GSS))


if __name__ == "__main__":
    unittest.main()
