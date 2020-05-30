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

import unittest

import numpy as np

import fwdpy11


class testGeneticValueIsFitness(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.g = fwdpy11.GeneticValueIsFitness(1)

    def testPickle(self):
        import pickle

        p = pickle.dumps(self.g, -1)
        up = pickle.loads(p)
        self.assertEqual(type(up), type(self.g))
        self.assertEqual(self.g.shape, up.shape)


class testGSS(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.VS = 1.0
        self.opt = 0.0
        self.g = fwdpy11.GSS(VS=self.VS, optimum=self.opt)

    def testPickle(self):
        import pickle

        p = pickle.dumps(self.g, -1)
        up = pickle.loads(p)
        self.assertEqual(up.VS, self.VS)
        self.assertEqual(up.optimum, self.opt)

    def test_mapping(self):
        self.assertEqual(self.g.maps_to_fitness, False)
        self.assertEqual(self.g.maps_to_trait_value, True)


class testGSSmo(unittest.TestCase):
    @classmethod
    def setUp(self):
        Opt = fwdpy11.Optimum
        self.optima = [
            Opt(when=0, optimum=0.0, VS=1.0),
            Opt(when=100, optimum=1.0, VS=1.0),
        ]
        self.g = fwdpy11.GSSmo(self.optima)

    def testPickle(self):
        import pickle

        p = pickle.dumps(self.g)
        up = pickle.loads(p)
        # self.assertEqual(up.VS, 1.0)
        # self.assertEqual(up.optimum, 0.0)
        for i, j in zip(up.optima, self.optima):
            self.assertEqual(i.when, j.when)
            self.assertEqual(i.optimum, j.optimum)
            self.assertEqual(i.VS, j.VS)

    def test_mapping(self):
        self.assertEqual(self.g.maps_to_fitness, False)
        self.assertEqual(self.g.maps_to_trait_value, True)


class testMultivariateGSSmo(unittest.TestCase):
    @classmethod
    def setUp(self):
        optima = np.array([0.0, 0.0, 1.0, 1.0]).reshape(2, 2)
        timepoints = [0, 100]
        po = []
        for i, t in enumerate(timepoints):
            po.append(fwdpy11.PleiotropicOptima(when=t, optima=optima[i, :], VS=1.0))
        self.mvgssmo = fwdpy11.MultivariateGSSmo(po)

    def test_pickle(self):
        import pickle

        p = pickle.dumps(self.mvgssmo)
        up = pickle.loads(p)
        self.assertEqual(self.mvgssmo, up)

    def test_mapping(self):
        self.assertEqual(self.mvgssmo.maps_to_fitness, False)
        self.assertEqual(self.mvgssmo.maps_to_trait_value, True)


if __name__ == "__main__":
    unittest.main()
