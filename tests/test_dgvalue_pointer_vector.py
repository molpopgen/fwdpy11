import unittest

import numpy as np

import fwdpy11


class TestConstruction(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.a = fwdpy11.Additive(2)
        self.list = [fwdpy11.Additive(2), fwdpy11.Additive(2)]

    def test_construct_from_object(self):
        from fwdpy11._fwdpy11 import _dgvalue_pointer_vector

        _dgvalue_pointer_vector(self.a)

    def test_construct_from_list(self):
        from fwdpy11._fwdpy11 import _dgvalue_pointer_vector

        _dgvalue_pointer_vector(self.list)


class TestBadInput(unittest.TestCase):
    def test_construct_from_empty_list(self):
        from fwdpy11._fwdpy11 import _dgvalue_pointer_vector

        with self.assertRaises(ValueError):
            _dgvalue_pointer_vector([])

    def test_invalid_dimensionality(self):
        a = fwdpy11.Additive(2.0)
        mvgss = fwdpy11.MultivariateGSS(np.zeros(2), 10.0)
        mv = fwdpy11.StrictAdditiveMultivariateEffects(2, 0, mvgss)
        from fwdpy11._fwdpy11 import _dgvalue_pointer_vector

        with self.assertRaises(ValueError):
            _dgvalue_pointer_vector([a, mv])

    def test_invalid_types(self):
        from fwdpy11._fwdpy11 import _dgvalue_pointer_vector

        with self.assertRaises(RuntimeError):
            _dgvalue_pointer_vector([1, 2])


if __name__ == "__main__":
    unittest.main()
