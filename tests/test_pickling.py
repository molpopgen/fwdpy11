#
# Copyright (C) 2017 Kevin Thornton <krthornt@uci.edu>
#
# This file is part of fwdpy11.
#
# fwdpy11 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# fwdpy11 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with fwdpy11.  If not, see <http://www.gnu.org/licenses/>.
#
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
        x = pickle.dumps(self.pop.mutations)
        p = pickle.loads(x)
        for i, j in zip(p, self.pop.mutations):
            self.assertEqual(i, j)

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

    def testPicklePopCpp(self):
        import pickle
        p = pickling_cpp.general_pickler(self.pop)
        pp = pickle.loads(p)
        self.assertEqual(pp, self.pop)


class testPickleMlocusPop(unittest.TestCase):
    @classmethod
    def setUp(self):
        from quick_pops import quick_mlocus_qtrait
        self.pop = quick_mlocus_qtrait(N=500, simlen=10)

    def testPickleDiploids(self):
        import pickle
        x = pickle.dumps(self.pop.diploids)
        p = pickle.loads(x)
        for i, j in zip(p, self.pop.diploids):
            self.assertEqual(i, j)

    def testPicklePopPy(self):
        import pickle
        p = pickle.dumps(self.pop, -1)
        pp = pickle.loads(p)
        self.assertEqual(pp, self.pop)

    def testPicklePopCpp(self):
        import pickle
        p = pickling_cpp.general_pickler(self.pop)
        pp = pickle.loads(p)
        self.assertEqual(pp, self.pop)


# class testSlocusPopGeneralMut(unittest.TestCase):
#     @classmethod
#     def setUp(self):
#         import fwdpy11
#         self.pop = fwdpy11.SlocusPopGeneralMutVec(100)
# 
#     def testPicklePopPy(self):
#         import pickle
#         p = pickle.dumps(self.pop, -1)
#         pp = pickle.loads(p)
#         self.assertEqual(pp, self.pop)
# 
#     def testPicklePopCpp(self):
#         import pickle
#         p = pickling_cpp.general_pickler(self.pop)
#         pp = pickle.loads(p)
#         self.assertEqual(pp, self.pop)
# 
# 
# class testVectorGeneralMutVec(unittest.TestCase):
#     """
#     As of 0.1.4, SlocusPopGeneralMutVec is
#     not used by any sims, and so there
#     is no evolve function to make a pop. Therefore,
#     we test that vectors of the mutation type
#     are pickleable.  We test that by making one.
#     """
#     @classmethod
#     def setUp(self):
#         import fwdpy11
#         s = fwdpy11.VecDouble([0.1, -0.1])
#         h = fwdpy11.VecDouble([0., 1.])
#         self.mutations = [fwdpy11.GeneralMutVec((s, h, 0.1, 3, 0))]
# 
#     def testPickleVector(self):
#         import pickle
#         p = pickle.dumps(self.mutations)
#         pp = pickle.loads(p)
#         for i, j in zip(pp, self.mutations):
#             self.assertEqual(i, j)


if __name__ == "__main__":
    unittest.main()
