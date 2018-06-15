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
import fwdpy11.internal
import numpy as np


class testMutationRegionWeights(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.spop = fwdpy11.SlocusPop(100)
        self.mpop = fwdpy11.MlocusPop(100, [(0, 100), (300, 1000)])
        self.mutrate_n = 1e-5
        self.mutrate_s = 1e-3
        self.pneutral = self.mutrate_n/(self.mutrate_n+self.mutrate_s)
        self.rng = fwdpy11.GSLrng(42)

    def testUniformWeightsSlocus(self):
        """
        Within neutral, selected, all regions are equally-likely to be picked.
        """
        nregions = [fwdpy11.Region(0, 1, 1), fwdpy11.Region(1, 2, 1)]
        sregions = [fwdpy11.ExpS(0, 1, 1, -0.1),
                    fwdpy11.ExpS(1, 2, 1, -0.2),
                    fwdpy11.ConstantS(0, 2, 1, -0.5, 1.0, False)]
        mm = fwdpy11.internal.makeMutationRegions(
            self.rng, self.spop, nregions, sregions, self.pneutral)
        w = np.array(mm.weights)
        pn = w[:2].sum()/w.sum()
        ps = w[2:].sum()/w.sum()
        self.assertEqual(pn, self.pneutral)
        self.assertEqual(ps, 1.0 - self.pneutral)

    def testWeightByRegionSizeSlocus(self):
        """
        Mutation rates are "per base pair", and the total weight that gets assigned 
        to a region is (end-beg)*mutation rate.  Thus Pr(region|mutation) is
        according to a region's length.
        """
        nregion_starts = [0, 100, 5000]
        nregion_stops = [100, 5000, 1e6]  # b/c intervals are [start,stop)
        nregion_lengths = np.array([b-a for a, b in zip(nregion_starts, nregion_stops)],
                                   dtype=np.float64)

        nregions = [fwdpy11.Region(start, stop, self.mutrate_n)
                    for start, stop in zip(nregion_starts, nregion_stops)]

        sregions = [fwdpy11.ExpS(start, stop, self.mutrate_s, -0.1)
                    for start, stop in zip(nregion_starts, nregion_stops)]

        mm = fwdpy11.internal.makeMutationRegions(
            self.rng, self.spop, nregions, sregions, self.pneutral)

        w = np.array(mm.weights)

        pn = w[:len(nregions)].sum()/w.sum()
        ps = w[len(nregions):].sum()/w.sum()

        # We assert almost equal b/c a+b != b+a
        self.assertAlmostEqual(pn, self.pneutral)
        self.assertAlmostEqual(ps, 1.0 - self.pneutral)

        # Check Pr(region x | mutation is neutral)
        nw = w[:len(nregions)]
        for i in range(len(nw)):
            pi = nw[i]/nw.sum()
            ei = nregion_lengths[i]/nregion_lengths.sum()
            self.assertAlmostEqual(pi,ei)

        # Check Pr(region x | mutation is selected)
        sw = w[len(nregions):]
        for i in range(len(nw)):
            pi = sw[i]/sw.sum()
            ei = nregion_lengths[i]/nregion_lengths.sum()
            self.assertAlmostEqual(pi,ei)


if __name__ == "__main__":
    unittest.main()
