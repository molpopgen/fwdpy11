import unittest
import msprime
import fwdpy11.ts_from_msprime
import numpy as np


class testConversion(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.ts = msprime.simulate(10, recombination_rate=0.025, Ne=1000)

    def testGetTables(self):
        n, e = fwdpy11.ts_from_msprime.convert_tables(self.ts, False)
        na = np.array(n, copy=False)
        self.assertTrue(np.all(na['time'][:10] == 0))

    def testGetTablesDiscretizeTime(self):
        n, e = fwdpy11.ts_from_msprime.convert_tables(self.ts, True)
        na = np.array(n, copy=False)
        tzero = np.where(na['time'] == 0.0)
        self.assertTrue(len(tzero[0]) == 10)


if __name__ == "__main__":
    unittest.main()
