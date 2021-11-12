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


class TestDiploidPopulation(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.pop = fwdpy11.DiploidPopulation(1000, 1.0)

    def test_N(self):
        self.assertEqual(self.pop.N, 1000)

    def test_founder_genome_counts(self):
        self.assertEqual(1, len(self.pop.haploid_genomes))
        self.assertEqual(2 * self.pop.N, self.pop.haploid_genomes[0].n)

    def test_generation(self):
        self.assertEqual(self.pop.generation, 0)

    def test_fitnesses(self):
        # All fitnesses should be 1
        self.assertEqual(
            sum([i.w for i in self.pop.diploid_metadata]), float(self.pop.N)
        )

    def test_labels(self):
        # All diploids should be labeled 0 to pop.N-1
        self.assertEqual(
            [i.label for i in self.pop.diploid_metadata], [i for i in range(self.pop.N)]
        )

    def test_genetic_values(self):
        self.assertEqual(sum([i.g for i in self.pop.diploid_metadata]), 0.0)

    def test_e_values(self):
        self.assertEqual(sum([i.e for i in self.pop.diploid_metadata]), 0.0)


class TestDiploidPopulationExceptions(unittest.TestCase):
    def testNzero(self):
        with self.assertRaises(ValueError):
            fwdpy11.DiploidPopulation(0, 1.0)


# NOTE: this tests an internal function (pop._record_ancient_samples),
# which is NOT the canonical user-facing API for doing this!
def test_ancient_sample_recording():
    pop = fwdpy11.DiploidPopulation(1000, 1.0)
    pop._record_ancient_samples([i for i in range(pop.N)])
    md = np.array(pop.diploid_metadata, copy=False)
    amd = np.array(pop.ancient_sample_metadata, copy=False)
    assert np.array_equal(md, amd)

    pop = fwdpy11.DiploidPopulation(1000, 1.0)
    pop._record_ancient_samples([i for i in range(0, pop.N, 2)])
    md = np.array(pop.diploid_metadata, copy=False)
    amd = np.array(pop.ancient_sample_metadata, copy=False)
    assert len(amd) == pop.N // 2
    assert np.array_equal(md[::2], amd)


if __name__ == "__main__":
    unittest.main()
