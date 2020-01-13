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
import unittest
import fwdpy11
import numpy as np


class TestSimpleMovesAndCopies(unittest.TestCase):
    """
    Basically doing the same tests
    as test_discrete_degraphy.TestDiscreteDemography
    """
    @classmethod
    def setUp(self):
        self.rng = fwdpy11.GSLrng(42)
        self.pop = fwdpy11.DiploidPopulation(100, 1.)
        self.pdict = {'nregions': [],
                      'sregions': [],
                      'recregions': [],
                      'rates': (0, 0, 0),
                      'demography': None,
                      'simlen': 1,
                      'gvalue': fwdpy11.Additive(2.)
                      }

    def test_simple_moves_from_single_deme(self):
        d = fwdpy11.DiscreteDemography(
            [fwdpy11.move_individuals(0, 0, 1, 0.5)])
        self.pdict['demography'] = d
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)
        md = np.array(self.pop.diploid_metadata, copy=False)
        deme_counts = np.unique(md['deme'], return_counts=True)
        self.assertEqual(len(deme_counts[0]), 2)
        for i in range(len(deme_counts[0])):
            self.assertEqual(deme_counts[1][i], 50)


if __name__ == "__main__":
    unittest.main()
