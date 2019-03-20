# Copyright (C) 2017 Kevin Thornton <krthornt@uci.edu>
#
# This file is part of fwdpy11.
#
# fwdpy11 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your meanion) any later version.
#
# fwdpy11 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with fwdpy11.  If not, see <http://www.gnu.org/licenses/>.

# Main goal is to test pickling

import fwdpy11
import unittest
import numpy as np


class testGeneticValueIsFitness(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.g = fwdpy11.GeneticValueIsFitness()

    def testPickle(self):
        import pickle
        p = pickle.dumps(self.g, -1)
        up = pickle.loads(p)
        self.assertEqual(type(up), type(self.g))


class testGSS(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.VS = 1.0
        self.opt = 0.0
        self.g = fwdpy11.GSS(VS=self.VS, opt=self.opt)

    def testPickle(self):
        import pickle
        p = pickle.dumps(self.g, -1)
        up = pickle.loads(p)
        self.assertEqual(up.VS, self.VS)
        self.assertEqual(up.opt, self.opt)


class testGSSmo(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.optima = [(0, 0.0, 1.0), (100, 1.0, 1.0)]
        self.g = fwdpy11.GSSmo(self.optima)

    def testPickle(self):
        import pickle
        p = pickle.dumps(self.g)
        up = pickle.loads(p)
        self.assertEqual(up.VS, 1.0)
        self.assertEqual(up.opt, 0.0)
        self.assertEqual(up.optima, self.optima)


class testMultivariateGSSmo(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.optima = np.array([0., 0., 1., 1.]).reshape(2, 2)
        self.timepoints = [0, 100]
        self.mvgssmo = fwdpy11.MultivariateGSSmo(
            self.timepoints, self.optima, 1.0)

    def test_pickle(self):
        import pickle
        p = pickle.dumps(self.mvgssmo)
        up = pickle.loads(p)
        self.assertEqual(self.mvgssmo, up)


if __name__ == "__main__":
    unittest.main()
