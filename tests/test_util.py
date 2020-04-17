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

# Unit tests of fwdpy11

import unittest

import fwdpy11
from quick_pops import quick_neutral_slocus


class testAddMutations(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.pop = quick_neutral_slocus()
        self.rng = fwdpy11.GSLrng(42)


class test_ChangeEsizeDiploid(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.pop = quick_neutral_slocus()
        self.rng = fwdpy11.GSLrng(42)

    def test_change_esize(self):
        # Get extant mutations
        extant = [i for i in enumerate(self.pop.mcounts) if i[1] > 0]
        # Get gametes with this mutation
        g = [
            i[0]
            for i in enumerate(self.pop.haploid_genomes)
            if extant[0][0] in i[1].mutations and self.pop.haploid_genomes[i[0]].n > 0
        ]
        self.assertEqual(sum([self.pop.haploid_genomes[i].n for i in g]), extant[0][1])
        fwdpy11.change_effect_size(self.pop, extant[0][0], -0.1)
        self.assertEqual(self.pop.mutations[extant[0][0]].s, -0.1)
        self.assertEqual(self.pop.mutations[extant[0][0]].neutral, False)
        g2 = [
            i[0]
            for i in enumerate(self.pop.haploid_genomes)
            if extant[0][0] in i[1].smutations and self.pop.haploid_genomes[i[0]].n > 0
        ]
        self.assertEqual(g, g2)

    def test_change_esize_nothing_happens(self):
        """
        Same test as above, but we change s from 0 to 0
        """
        # Get extant mutations
        extant = [i for i in enumerate(self.pop.mcounts) if i[1] > 0]
        # Get gametes with this mutation
        g = [
            i[0]
            for i in enumerate(self.pop.haploid_genomes)
            if extant[0][0] in i[1].mutations and self.pop.haploid_genomes[i[0]].n > 0
        ]
        self.assertEqual(sum([self.pop.haploid_genomes[i].n for i in g]), extant[0][1])
        fwdpy11.change_effect_size(self.pop, extant[0][0], 0.0)
        self.assertEqual(self.pop.mutations[extant[0][0]].s, 0.0)
        self.assertEqual(self.pop.mutations[extant[0][0]].neutral, True)
        g2 = [
            i[0]
            for i in enumerate(self.pop.haploid_genomes)
            if extant[0][0] in i[1].mutations and self.pop.haploid_genomes[i[0]].n > 0
        ]
        self.assertEqual(g, g2)

    def testVectorEffects(self):
        extant = [i for i in enumerate(self.pop.mcounts) if i[1] > 0]
        g = [
            i[0]
            for i in enumerate(self.pop.haploid_genomes)
            if extant[0][0] in i[1].mutations and self.pop.haploid_genomes[i[0]].n > 0
        ]
        fwdpy11.change_effect_size(self.pop, extant[0][0], 0.0, 0.0, [1.0], [-1.0])
        self.assertEqual(self.pop.mutations[extant[0][0]].neutral, False)
        self.assertEqual(self.pop.mutations[extant[0][0]].esizes[0], 1.0)
        self.assertEqual(self.pop.mutations[extant[0][0]].heffects[0], -1.0)

        # Now, replace it w/a neutral mutation
        fwdpy11.change_effect_size(self.pop, extant[0][0], 0.0, 0.0, [0.0], [-1.0])
        self.assertEqual(self.pop.mutations[extant[0][0]].neutral, True)
        self.assertEqual(self.pop.mutations[extant[0][0]].esizes[0], 0.0)
        self.assertEqual(self.pop.mutations[extant[0][0]].heffects[0], -1.0)


if __name__ == "__main__":
    unittest.main()
