import unittest
import fwdpy11
import numpy as np


class testConversion(unittest.TestCase):
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
        n = md['nodes'].flatten()
        self.assertTrue(np.array_equal(
            n, np.arange(2*pop.N, dtype=n.dtype)))


if __name__ == "__main__":
    unittest.main()
