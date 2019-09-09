import unittest
import numpy_array_interface


class TestNumpyArrayInterface(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.v = numpy_array_interface.VecBackedArray()

    def test_readwrite(self):
        a = self.v.x_readwrite()
        self.assertEqual(a.flags.owndata, False)
        self.assertEqual(a.flags.writeable, True)

    def test_readonly(self):
        a = self.v.x_readonly()
        self.assertEqual(a.flags.owndata, False)
        self.assertEqual(a.flags.writeable, False)


if __name__ == "__main__":
    unittest.main()
