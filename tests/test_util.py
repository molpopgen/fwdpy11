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

# Unit tests of fwdpy11.util

import unittest
from quick_pops import quick_neutral_slocus, quick_mlocus_qtrait
import fwdpy11
import fwdpy11.util


class testAddMutations(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.pop = quick_neutral_slocus()
        self.rng = fwdpy11.GSLrng(42)

    def test_add_selected(self):
        ncopies = 103
        x = fwdpy11.util.add_mutation(
            self.rng, self.pop, ncopies, (0.77, -0.1, 0.11))
        self.assertTrue(x < len(self.pop.mutations))
        self.assertTrue(self.pop.mcounts[x] == ncopies)

    def test_mutate_at_existing_position(self):
        ncopies = 213
        extant_pos = [self.pop.mutations[i].pos
                      for i in range(len(self.pop.mcounts))
                      if self.pop.mcounts[i] > 0]
        with self.assertRaises(ValueError):
            x = fwdpy11.util.add_mutation(
                self.rng, self.pop, ncopies, (extant_pos[0], -0.1, 0.11))
            x


class test_ChangeEsizeSlocus(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.pop = quick_neutral_slocus()
        self.rng = fwdpy11.GSLrng(42)

    def test_change_esize(self):
        self.setUpClass()
        # Get extant mutations
        extant = [i for i in enumerate(self.pop.mcounts) if i[1] > 0]
        # Get gametes with this mutation
        g = [i[0] for i in enumerate(
            self.pop.gametes) if extant[0][0] in i[1].mutations and
            self.pop.gametes[i[0]].n > 0]
        self.assertEqual(sum([self.pop.gametes[i].n for i in g]), extant[0][1])
        fwdpy11.util.change_effect_size(self.pop, extant[0][0], -0.1)
        self.assertEqual(self.pop.mutations[extant[0][0]].s, -0.1)
        self.assertEqual(self.pop.mutations[extant[0][0]].neutral, False)
        g2 = [i[0] for i in enumerate(self.pop.gametes) if
              extant[0][0] in i[1].smutations and
              self.pop.gametes[i[0]].n > 0]
        self.assertEqual(g, g2)

    def test_change_esize_nothing_happens(self):
        """
        Same test as above, but we change s from 0 to 0
        """
        self.setUpClass()
        # Get extant mutations
        extant = [i for i in enumerate(self.pop.mcounts) if i[1] > 0]
        # Get gametes with this mutation
        g = [i[0] for i in enumerate(
            self.pop.gametes) if extant[0][0] in i[1].mutations and
            self.pop.gametes[i[0]].n > 0]
        self.assertEqual(sum([self.pop.gametes[i].n for i in g]), extant[0][1])
        fwdpy11.util.change_effect_size(self.pop, extant[0][0], 0.0)
        self.assertEqual(self.pop.mutations[extant[0][0]].s, 0.0)
        self.assertEqual(self.pop.mutations[extant[0][0]].neutral, True)
        g2 = [i[0] for i in enumerate(self.pop.gametes) if
              extant[0][0] in i[1].mutations and
              self.pop.gametes[i[0]].n > 0]
        self.assertEqual(g, g2)


class test_ChangeEsizeMlocus(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.pop = quick_mlocus_qtrait()

    def test_change_esize(self):
        extant = [i for i in enumerate(self.pop.mcounts) if i[1] > 0]
        # Get gametes with this mutation
        g = [i[0] for i in enumerate(
            self.pop.gametes) if extant[0][0] in i[1].mutations and
            self.pop.gametes[i[0]].n > 0]
        self.assertEqual(sum([self.pop.gametes[i].n for i in g]), extant[0][1])
        fwdpy11.util.change_effect_size(self.pop, extant[0][0], -0.1)
        self.assertEqual(self.pop.mutations[extant[0][0]].s, -0.1)
        self.assertEqual(self.pop.mutations[extant[0][0]].neutral, False)
        g2 = [i[0] for i in enumerate(self.pop.gametes) if
              extant[0][0] in i[1].smutations and
              self.pop.gametes[i[0]].n > 0]
        self.assertEqual(g, g2)


if __name__ == "__main__":
    unittest.main()
