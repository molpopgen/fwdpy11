import fwdpy11
import unittest


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


if __name__ == "__main__":
    unittest.main()
