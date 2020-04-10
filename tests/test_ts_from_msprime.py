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

import numpy as np

import fwdpy11

# NOTE: msprime is imported on a per-test basis so as not
# to confict w/other tests of the GSL error handling machinery


class TestConversion(unittest.TestCase):
    @classmethod
    def setUp(self):
        import msprime

        self.ts = msprime.simulate(10, recombination_rate=0.025, Ne=1000)

    # def testGetTablesDiscretizeTime(self):
    #     n, e, ntips, l = fwdpy11.ts_from_msprime._convert_tables(self.ts)
    #     self.assertEqual(ntips, 10)
    #     self.assertEqual(l, 1.0)
    #     na = np.array(n, copy=False)
    #     tzero = np.where(na['time'] == 0.0)
    #     self.assertTrue(len(tzero[0]) == 10)

    def testCreateDiploidPopulation(self):
        pop = fwdpy11.DiploidPopulation.create_from_tskit(self.ts)
        self.assertEqual(pop.N, 5)
        self.assertEqual(pop.tables.genome_length, 1.0)
        md = np.array(pop.diploid_metadata, copy=False)
        n = md["nodes"].flatten()
        self.assertTrue(np.array_equal(n, np.arange(2 * pop.N, dtype=n.dtype)))


class TestConversionFromMultipleDemes(unittest.TestCase):
    def test_deme_field_of_metadata(self):
        import msprime

        nodes_per_deme = 500
        Ne = 3 * nodes_per_deme // 2
        Nr = 50.0
        config = [
            msprime.PopulationConfiguration(nodes_per_deme),
            msprime.PopulationConfiguration(nodes_per_deme),
            msprime.PopulationConfiguration(nodes_per_deme),
        ]

        events = [
            msprime.MassMigration(1 * Ne, 1, 0, 1.0),
            msprime.MassMigration(1.5 * Ne, 2, 0, 1.0),
        ]

        ts = msprime.simulate(
            population_configurations=config,
            Ne=Ne,
            random_seed=98765,
            recombination_rate=Nr / Ne,
            demographic_events=events,
        )
        pop = fwdpy11.DiploidPopulation.create_from_tskit(ts)
        self.assertEqual(pop.N, 750)
        alive_nodes = pop.alive_nodes
        self.assertEqual(len(alive_nodes), 2 * pop.N)
        for i in range(3):
            for j in alive_nodes[i * nodes_per_deme : (i + 1) * nodes_per_deme]:
                self.assertEqual(pop.tables.nodes[j].deme, i)
            k, el = i * nodes_per_deme // 2, (i + 1) * nodes_per_deme // 2
            for j in pop.diploid_metadata[k:el]:
                self.assertEqual(j.deme, i)


if __name__ == "__main__":
    unittest.main()
