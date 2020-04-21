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

import pickle
import unittest

import fwdpy11


class testSingleEffectMutation(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.m = fwdpy11.Mutation(1.0, -1.0, 0.25, 0, 13)

    def testConstruct(self):
        self.assertEqual(self.m.pos, 1.0)
        self.assertEqual(self.m.s, -1.0)
        self.assertEqual(self.m.h, 0.25)
        self.assertEqual(self.m.g, 0)
        self.assertEqual(self.m.label, 13)
        self.assertEqual(len(self.m.esizes), 0)
        self.assertEqual(len(self.m.heffects), 0)

    def testPickle(self):
        p = pickle.dumps(self.m)
        up = pickle.loads(p)
        self.assertEqual(up, self.m)


class testMultiEffectMutation(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.m = fwdpy11.Mutation(1.0, -1.0, 0.25, 0, [1.0, 2.0], [-1.0, 1.0 / 3.0], 13)

    def testConstruct(self):
        self.assertEqual(self.m.pos, 1.0)
        self.assertEqual(self.m.s, -1.0)
        self.assertEqual(self.m.h, 0.25)
        self.assertEqual(self.m.g, 0)
        self.assertEqual(self.m.label, 13)
        self.assertEqual(len(self.m.esizes), 2)
        self.assertEqual(len(self.m.heffects), 2)

    def testPickle(self):
        p = pickle.dumps(self.m)
        up = pickle.loads(p)
        self.assertEqual(up, self.m)


if __name__ == "__main__":
    unittest.main()
