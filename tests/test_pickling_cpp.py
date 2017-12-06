import cppimport
cppimport.force_rebuild()
cppimport.set_quiet(False)
pickling_cpp = cppimport.imp("pickling_cpp")
import unittest


class testPickeMutation(unittest.TestCase):
    @classmethod
    def setUp(self):
        import fwdpy11
        self.m = fwdpy11.Mutation(0.1, -0.2, 1.0, 10, 3)

    def testUnpickle(self):
        import pickle
        o = pickling_cpp.pickle_mutation(self.m)
        m = pickle.loads(o)
        self.assertEqual(m, self.m)


if __name__ == "__main__":
    unittest.main()
