import unittest

import pickling_composed_classes as PC


class testPickleComposedClass(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.D = PC.Derived()
        self.C = PC.Composed(2.0, self.D)

    def testPickle(self):
        import pickle

        p = pickle.dumps(self.C, -1)
        up = pickle.loads(p)
        self.assertEqual(type(up), type(self.C))
        self.assertEqual(type(up.b), type(self.C.b))


if __name__ == "__main__":
    unittest.main()
