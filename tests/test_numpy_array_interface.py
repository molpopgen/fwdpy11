import unittest

import numpy_array_interface


class TestNumpyArrayInterface1D(unittest.TestCase):
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

    def test_readwrite_via_capsule(self):
        a = self.v.x_readwrite_via_capsule()
        self.assertEqual(a.flags.owndata, False)
        self.assertEqual(a.flags.writeable, True)
        b = self.v.x_readwrite()
        self.assertEqual(b.shape[0], 0)


class TestNumpyArrayInterface2D(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.v = numpy_array_interface.VecBacked2DArray(7, 3)

    def test_shape(self):
        a = self.v.x_readwrite()
        self.assertEqual(a.shape, self.v.shape)

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
