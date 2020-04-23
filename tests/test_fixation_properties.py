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
# Test behavior of fixations in simulations.
# As of 0.1.3a1, fixations are to be sorted by position
# for all types of simulations, allowing for easier lookup.
#
# In 0.2.0, this unit test was refactored to skip simulation.
# Instead, we use pybind11 + cppimport to expose the underlying
# C++ API to Python, create pops manually, and test the API.

import os
import unittest

import fixation_properties as fp
import fwdpy11


@unittest.skip("Test uses removed API and needs updating")
class testFixationsAreSortedDiploidPopulation(unittest.TestCase):
    @classmethod
    def setUp(self):
        mutations = fwdpy11.MutationVector()
        fixations = fwdpy11.MutationVector()
        gametes = fwdpy11.HaploidGenomeVector()
        diploids = fwdpy11.DiploidVector()
        mutations.append(fwdpy11.Mutation(0.1, -0.01, 1.0, 0, 0))
        mutations.append(fwdpy11.Mutation(0.2, 0.0, 1.0, 0, 0))
        gametes.append(fwdpy11.HaploidGenome((4, [1], [0])))
        diploids.append(fwdpy11.DiploidGenotype(0, 0))
        diploids.append(fwdpy11.DiploidGenotype(0, 0))
        self.pop = fwdpy11.DiploidPopulation(diploids, gametes, mutations)

    def testSetup(self):
        self.assertEqual(len(self.pop.mutations), 2)
        self.assertEqual(len(self.pop.haploid_genomes), 1)
        self.assertEqual(len(self.pop.diploids), 2)
        for g in self.pop.haploid_genomes:
            self.assertEqual(len(g.mutations), 1)
            self.assertEqual(len(g.smutations), 1)
        for mc in self.pop.mcounts:
            self.assertEqual(mc, 4)

    def testWithPruneSelected(self):
        fp.update_mutations(self.pop, 1, 2 * len(self.pop.diploids), True)
        self.assertEqual(len(self.pop.fixations), 2)
        self.assertEqual(len(self.pop.fixation_times), 2)
        self.assertEqual(len([i for i in self.pop.mcounts if i == 0]), 2)
        self.assertTrue(sorted(self.pop.fixations, key=lambda x: x.pos))
        for mc in self.pop.mcounts:
            self.assertEqual(mc, 0)

    def testWithoutPruneSelected(self):
        fp.update_mutations(self.pop, 1, 2 * len(self.pop.diploids), False)
        self.assertEqual(len(self.pop.fixations), 2)
        self.assertEqual(len(self.pop.fixation_times), 2)
        self.assertEqual(len([i for i in self.pop.mcounts if i == 0]), 1)
        self.assertTrue(sorted(self.pop.fixations, key=lambda x: x.pos))
        for m, mc in zip(self.pop.mutations, self.pop.mcounts):
            if m.neutral is True:
                self.assertEqual(mc, 0)
            else:
                self.assertEqual(mc, 4)


class testFixationPreservation(unittest.TestCase):
    @classmethod
    def setUp(self):
        import numpy as np

        N = 1000
        demography = np.array([N] * 10 * N, dtype=np.uint32)
        rho = 1.0
        r = rho / (4 * N)

        a = fwdpy11.Multiplicative(2.0)
        self.p = {
            "nregions": [],
            "sregions": [fwdpy11.ExpS(0, 1, 1, 0.01)],
            "recregions": [fwdpy11.Region(0, 1, 1)],
            "rates": (0.0, 0.00005, r),
            "gvalue": a,
            "popsizes": demography,
            "simlen": len(demography),
        }
        self.pop = fwdpy11.DiploidPopulation(N)
        self.rng = fwdpy11.GSLrng(101 * 45 * 110 * 210)

    def testPopGenSimWithoutPruning(self):
        import fwdpy11
        import numpy as np

        self.p["prune_selected"] = False
        params = fwdpy11.ModelParams(**self.p)
        fwdpy11.evolve_genomes(self.rng, self.pop, params)
        assert len(self.pop.fixations) > 0, "Test is meaningless without fixations"
        mc = np.array(self.pop.mcounts)
        self.assertEqual(
            len(np.where(mc == 2 * self.pop.N)[0]), len(self.pop.fixations)
        )

    def testPopGenSimWithPruning(self):
        import fwdpy11
        import numpy as np

        self.p["prune_selected"] = True
        params = fwdpy11.ModelParams(**self.p)
        fwdpy11.evolve_genomes(self.rng, self.pop, params)
        assert len(self.pop.fixations) > 0, "Test is meaningless without fixations"
        mc = np.array(self.pop.mcounts)
        self.assertEqual(len(np.where(mc == 2 * self.pop.N)[0]), 0)


if __name__ == "__main__":
    unittest.main()
