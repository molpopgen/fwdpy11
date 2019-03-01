import unittest
import fwdpy11.ts_from_tskit
import numpy as np
import fwdpy11


class testConversion(unittest.TestCase):
    @classmethod
    def setUp(self):
        import msprime
        self.ts = msprime.simulate(10, recombination_rate=0.025, Ne=1000)

    def testGetTables(self):
        n, e, ntips, l = fwdpy11.ts_from_tskit._convert_tables(self.ts)
        self.assertEqual(ntips, 10)
        self.assertEqual(l, 1.0)
        na = np.array(n, copy=False)
        self.assertTrue(np.all(na['time'][:10] == 0))

    # def testGetTablesDiscretizeTime(self):
    #     n, e, ntips, l = fwdpy11.ts_from_msprime._convert_tables(self.ts)
    #     self.assertEqual(ntips, 10)
    #     self.assertEqual(l, 1.0)
    #     na = np.array(n, copy=False)
    #     tzero = np.where(na['time'] == 0.0)
    #     self.assertTrue(len(tzero[0]) == 10)

    def testCreateSlocusPop(self):
        pop = fwdpy11.SlocusPop.create_from_tskit(self.ts)
        self.assertEqual(pop.N, 5)
        self.assertEqual(pop.tables.genome_length(), 1.0)


if __name__ == "__main__":
    unittest.main()
