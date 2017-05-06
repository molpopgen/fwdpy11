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

#Unit tests of fwdpy11.util

import unittest
from quick_pops import quick_neutral_slocus
import fwdpy11
import fwdpy11.util

class testAddMutations(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.pop=quick_neutral_slocus()
        self.rng=fwdpy11.GSLrng(42)
    def test_add_selected(self):
        ncopies = 103
        x=fwdpy11.util.add_mutation(self.rng,self.pop,ncopies,(0.77,-0.1,0.11))
        self.assertTrue(x < len(self.pop.mutations))
        self.assertTrue(self.pop.mcounts[x] == ncopies)
    def test_mutate_at_existing_position(self):
        ncopies = 213
        extant_pos=[self.pop.mutations[i].pos for i in range(len(self.pop.mcounts)) if self.pop.mcounts[i]>0]
        with self.assertRaises(ValueError):
            x=fwdpy11.util.add_mutation(self.rng,self.pop,ncopies,(extant_pos[0],-0.1,0.11))

if __name__ == "__main__":
    unittest.main()

