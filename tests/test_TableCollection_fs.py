#
# Copyright (C) 2020 Kevin Thornton <krthornt@uci.edu>
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

import fwdpy11
import unittest
import numpy as np


class TestSingleDemeCase(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        import msprime
        Ne = 1000
        Nr = 100.0
        ts = msprime.simulate(2*Ne, Ne=Ne, recombination_rate=Nr/Ne)
        self.pop = fwdpy11.DiploidPopulation.create_from_tskit(ts)
        rng = fwdpy11.GSLrng(12343)
        fwdpy11.infinite_sites(rng, self.pop, Nr/Ne)

    def test_comparison_to_genotype_matrix(self):
        dm = fwdpy11.data_matrix_from_tables(self.pop.tables,
                                             self.pop.alive_nodes, True, False)
        gm = np.array(dm.neutral)
        gm_rc = np.sum(gm, axis=1)
        gm_uc = np.unique(gm_rc, return_counts=True)
        gm_fs = np.zeros(2*self.pop.N+1, dtype=np.int32)
        gm_fs[gm_uc[0]] = gm_uc[1]
        tc_fs = self.pop.tables.fs([self.pop.alive_nodes])
        self.assertTrue(np.array_equal(gm_fs[1:-1], tc_fs.data[1:-1]))

    def test_comparison_to_genotype_matrix_for_sample(self):
        A, B = 103, 210
        dm = fwdpy11.data_matrix_from_tables(self.pop.tables,
                                             self.pop.alive_nodes[A:B],
                                             True, False)
        gm = np.array(dm.neutral)
        gm_rc = np.sum(gm, axis=1)
        gm_uc = np.unique(gm_rc, return_counts=True)
        gm_fs = np.zeros(B-A+1, dtype=np.int32)
        gm_fs[gm_uc[0]] = gm_uc[1]
        tc_fs = self.pop.tables.fs([self.pop.alive_nodes[A:B]])
        self.assertTrue(np.array_equal(gm_fs[1:-1], tc_fs.data[1:-1]))

    def test_skipping_neutral_variants(self):
        tc_fs = self.pop.tables.fs(
            [self.pop.alive_nodes], include_neutral=False)
        self.assertEqual(tc_fs.sum(), 0)

    def test_empty_samples_list(self):
        with self.assertRaises(ValueError):
            self.pop.tables.fs([])

    def test_nodes_out_of_range(self):
        with self.assertRaises(ValueError):
            samples = np.array([len(self.pop.tables.nodes)], dtype=np.int32)
            self.pop.tables.fs([samples])


if __name__ == "__main__":
    unittest.main()
