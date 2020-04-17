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


class testNoNoise(unittest.TestCase):
    @classmethod
    def setUp(self):
        import fwdpy11

        self.n = fwdpy11.NoNoise()

    def testPickle(self):
        import pickle

        p = pickle.dumps(self.n, -1)
        up = pickle.loads(p)
        self.assertEqual(type(up), type(self.n))


class testNoNoisePackageAlias(unittest.TestCase):
    @classmethod
    def setUp(self):
        import fwdpy11 as fancyname

        self.n = fancyname.NoNoise()

    def testPickle(self):
        import pickle

        p = pickle.dumps(self.n, -1)
        up = pickle.loads(p)
        self.assertEqual(type(up), type(self.n))


class testGaussianNoise(unittest.TestCase):
    @classmethod
    def setUp(self):
        import fwdpy11

        self.sd = 1.0
        self.mean = 0.0
        self.n = fwdpy11.GaussianNoise(sd=self.sd, mean=self.mean)

    def testPickle(self):
        import pickle

        p = pickle.dumps(self.n, -1)
        up = pickle.loads(p)
        self.assertEqual(up.sd, self.sd)
        self.assertEqual(up.mean, self.mean)


if __name__ == "__main__":
    unittest.main()
