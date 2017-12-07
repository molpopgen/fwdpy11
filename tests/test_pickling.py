import cppimport
cppimport.force_rebuild()
cppimport.set_quiet(False)
pickling_cpp = cppimport.imp("pickling_cpp")
import unittest


class testPickleMutation(unittest.TestCase):
    @classmethod
    def setUp(self):
        import fwdpy11
        self.m = fwdpy11.Mutation(0.1, -0.2, 1.0, 10, 3)

    def testUnpickle(self):
        import pickle
        o = pickling_cpp.pickle_mutation(self.m)
        m = pickle.loads(o)
        self.assertEqual(m, self.m)

    def testUnpickleFromGeneralPickler(self):
        import pickle
        o = pickling_cpp.general_pickler(self.m)
        m = pickle.loads(o)
        self.assertEqual(m, self.m)


class testPickleSlocusPop(unittest.TestCase):
    @classmethod
    def setUp(self):
        from quick_pops import quick_neutral_slocus
        self.pop = quick_neutral_slocus()

    def testPickleMutations(self):
        import pickle
        p = [pickling_cpp.general_pickler(i) for i in self.pop.mutations]
        for i, j in zip(p, self.pop.mutations):
            m = pickle.loads(i)
            self.assertEqual(m, j)

    def testPickleGametes(self):
        import pickle
        p = [pickling_cpp.general_pickler(i) for i in self.pop.gametes]
        for i, j in zip(p, self.pop.gametes):
            m = pickle.loads(i)
            self.assertEqual(m, j)
        x = pickling_cpp.general_pickler(self.pop.gametes)
        p = pickle.loads(x)
        for i, j in zip(p, self.pop.gametes):
            self.assertEqual(i, j)

    def testPickleDiploidsPy(self):
        import pickle
        pd = pickle.dumps(self.pop.diploids)
        p = pickle.loads(pd)
        for i, j in zip(p, self.pop.diploids):
            self.assertEqual(i, j)

    def testPickleDiploidsCpp(self):
        import pickle
        p = [pickling_cpp.general_pickler(i) for i in self.pop.diploids]
        for i, j in zip(p, self.pop.diploids):
            m = pickle.loads(i)
            self.assertEqual(m, j)
        x = pickling_cpp.general_pickler(self.pop.diploids)
        p = pickle.loads(x)
        for i, j in zip(p, self.pop.diploids):
            self.assertEqual(i, j)


if __name__ == "__main__":
    unittest.main()
