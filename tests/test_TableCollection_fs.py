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

import unittest

import fwdpy11
import msprime
import numpy as np

# NOTE: these tests can all be improved in the future
# by having a conversion from msprime/tskit that lifts
# over mutations assuming that they are neutral.
# With that in place, we can compare to the msprime/tskit
# genotype matrix, which makes for a more independent test
# of correctness.


def fs_from_ndarray(gm):
    gm_rc = np.sum(gm, axis=1)
    gm_uc = np.unique(gm_rc, return_counts=True)
    gm_fs = np.zeros(gm.shape[1] + 1, dtype=np.int32)
    gm_fs[gm_uc[0]] += gm_uc[1]
    return gm_fs


class TestSingleDemeCase(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        Ne = 1000
        Nr = 100.0
        ts = msprime.simulate(
            2 * Ne, Ne=Ne, recombination_rate=Nr / Ne, random_seed=666
        )
        self.pop = fwdpy11.DiploidPopulation.create_from_tskit(ts)
        rng = fwdpy11.GSLrng(12343)
        fwdpy11.infinite_sites(rng, self.pop, Nr / Ne)

    def test_comparison_to_genotype_matrix(self):
        dm = fwdpy11.data_matrix_from_tables(
            self.pop.tables, self.pop.alive_nodes, True, False
        )
        gm = np.array(dm.neutral, copy=False)
        gm_fs = fs_from_ndarray(gm)
        tc_fs = self.pop.tables.fs([self.pop.alive_nodes])
        self.assertTrue(np.array_equal(gm_fs[1:-1], tc_fs.data[1:-1]))
        self.assertTrue(np.array_equal(gm_fs, tc_fs.data))

    def test_comparison_to_genotype_matrix_for_sample(self):
        A, B = 103, 210
        dm = fwdpy11.data_matrix_from_tables(
            self.pop.tables, self.pop.alive_nodes[A:B], True, False
        )
        gm = np.array(dm.neutral, copy=False)
        gm_fs = fs_from_ndarray(gm)
        tc_fs = self.pop.tables.fs([self.pop.alive_nodes[A:B]])
        self.assertTrue(np.array_equal(gm_fs[1:-1], tc_fs.data[1:-1]))

    def test_skipping_neutral_variants(self):
        tc_fs = self.pop.tables.fs([self.pop.alive_nodes], include_neutral=False)
        self.assertEqual(tc_fs.sum(), 0)

    def test_empty_samples_list(self):
        with self.assertRaises(ValueError):
            self.pop.tables.fs([])

    def test_nodes_out_of_range(self):
        with self.assertRaises(ValueError):
            samples = np.array([len(self.pop.tables.nodes)], dtype=np.int32)
            self.pop.tables.fs([samples])

    def test_two_overlapping_windows(self):
        tc_fs = self.pop.tables.fs([self.pop.alive_nodes])
        wfs = self.pop.tables.fs(
            [self.pop.alive_nodes], windows=[(0, 1.0 / 3.0), (1.0 / 3.0, 1.0)]
        )
        self.assertTrue(np.array_equal(wfs, tc_fs))

    def test_three_overlapping_windows(self):
        tc_fs = self.pop.tables.fs([self.pop.alive_nodes])
        wfs = self.pop.tables.fs(
            [self.pop.alive_nodes],
            windows=[(0, 1.0 / 3.0), (1.0 / 3.0, 2.0 / 3.0), (2.0 / 3.0, 1)],
        )
        self.assertTrue(np.array_equal(wfs, tc_fs))

    def test_random_number_of_windows(self):
        nw = np.random.randint(100, 300)
        lefts = np.arange(nw) / (nw - 1)
        windows = [(lefts[i], lefts[i + 1]) for i in range(len(lefts) - 1)]
        tc_fs = self.pop.tables.fs([self.pop.alive_nodes])
        wfs = self.pop.tables.fs([self.pop.alive_nodes], windows=windows)
        self.assertTrue(np.array_equal(wfs, tc_fs))

    def test_separated_windows(self):
        windows = [(0, 0.25), (0.66, 0.9)]
        dm = fwdpy11.data_matrix_from_tables(
            self.pop.tables, self.pop.alive_nodes, True, False
        )
        gm = np.array(dm.neutral, copy=False)
        gm_pos = np.array([self.pop.mutations[k].pos for k in dm.neutral_keys])
        gm_pos_in_windows = np.where(
            (gm_pos < 0.25) | ((gm_pos >= 0.66) & (gm_pos < 0.9))
        )[0]
        gm = gm[gm_pos_in_windows, :]
        gm_fs = fs_from_ndarray(gm)
        tc_fs = self.pop.tables.fs([self.pop.alive_nodes], windows=windows)
        self.assertTrue(np.array_equal(gm_fs, tc_fs))


class TestTwoDemeCase(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        Ne = 1000
        Nr = 100.0
        nodes_per_deme = 1000
        config = [
            msprime.PopulationConfiguration(nodes_per_deme),
            msprime.PopulationConfiguration(nodes_per_deme),
        ]

        events = [msprime.MassMigration(1 * Ne, 1, 0, 1.0)]
        ts = msprime.simulate(
            population_configurations=config,
            demographic_events=events,
            Ne=Ne,
            recombination_rate=Nr / Ne,
            random_seed=666,
        )
        self.pop = fwdpy11.DiploidPopulation.create_from_tskit(ts)
        rng = fwdpy11.GSLrng(12343)
        fwdpy11.infinite_sites(rng, self.pop, Nr / Ne)

    def test_marginal_deme_fs(self):
        a = self.pop.alive_nodes
        nodes = np.array(self.pop.tables.nodes, copy=False)
        d0 = a[np.where(nodes["deme"][a] == 0)[0]]
        d1 = a[np.where(nodes["deme"][a] == 1)[0]]

        for samples in (d0, d1):
            tc_fs = self.pop.tables.fs([samples])
            dm = fwdpy11.data_matrix_from_tables(self.pop.tables, samples, True, True)
            gm = np.array(dm.neutral, copy=False)
            gm_fs = fs_from_ndarray(gm)
            self.assertTrue(np.array_equal(gm_fs[1:-1], tc_fs.data[1:-1]))

    def test_joint_deme_fs(self):
        a = self.pop.alive_nodes
        nodes = np.array(self.pop.tables.nodes, copy=False)
        d0 = a[np.where(nodes["deme"][a] == 0)[0]]
        d1 = a[np.where(nodes["deme"][a] == 1)[0]]

        tc_fs = self.pop.tables.fs([d0, d1])
        tc_fs0 = tc_fs.sum(axis=1).todense()
        tc_fs1 = tc_fs.sum(axis=0).todense()
        for i, j in zip((tc_fs0, tc_fs1), (d0, d1)):
            dm = fwdpy11.data_matrix_from_tables(self.pop.tables, j, True, True)
            gm = np.array(dm.neutral, copy=False)
            gm_fs = fs_from_ndarray(gm)
            self.assertTrue(np.array_equal(gm_fs[1:-1], i.data[1:-1]))

    def test_joint_deme_fs_marginalize(self):
        a = self.pop.alive_nodes
        nodes = np.array(self.pop.tables.nodes, copy=False)
        d0 = a[np.where(nodes["deme"][a] == 0)[0]]
        d1 = a[np.where(nodes["deme"][a] == 1)[0]]
        tc_fs = self.pop.tables.fs([d0, d1])
        tc_fs_deme0 = tc_fs.sum(axis=1).todense()
        tc_fs_deme1 = tc_fs.sum(axis=0).todense()

        tc_fs2 = self.pop.tables.fs([d0, d1], marginalize=True)
        self.assertTrue(np.array_equal(tc_fs_deme0, tc_fs2[0]))
        self.assertTrue(np.array_equal(tc_fs_deme1, tc_fs2[1]))

    def test_marginalizing_to_three_samples(self):
        a = self.pop.alive_nodes
        sample_lists = [a[:10], a[20:30], a[200:210]]
        fs1 = self.pop.tables.fs([sample_lists[0]])
        fs2 = self.pop.tables.fs([sample_lists[1]])
        fs3 = self.pop.tables.fs([sample_lists[2]])

        mfs = self.pop.tables.fs(sample_lists, marginalize=True)

        for i, j in zip(range(3), [fs1, fs2, fs3]):
            self.assertTrue(np.ma.allequal(mfs[i], j))

    def test_marginalizing_to_three_samples_in_separate_windows(self):
        a = self.pop.alive_nodes
        sample_lists = [a[:10], a[20:30], a[200:210]]
        windows = [(0.1, 0.23), (0.723, 0.85)]
        fs1 = self.pop.tables.fs(
            [sample_lists[0]], windows=windows, separate_windows=True
        )
        fs2 = self.pop.tables.fs(
            [sample_lists[1]], windows=windows, separate_windows=True
        )
        fs3 = self.pop.tables.fs(
            [sample_lists[2]], windows=windows, separate_windows=True
        )

        mfs = self.pop.tables.fs(
            sample_lists, marginalize=True, windows=windows, separate_windows=True
        )

        for i, j in zip(mfs, range(len(windows))):
            self.assertTrue(np.ma.allequal(i[0], fs1[j]))
            self.assertTrue(np.ma.allequal(i[1], fs2[j]))
            self.assertTrue(np.ma.allequal(i[2], fs3[j]))


if __name__ == "__main__":
    unittest.main()
