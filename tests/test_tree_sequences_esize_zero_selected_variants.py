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

# This test is motivated by GitHub issue #432.
# When a DES/DFE returns an effect size of zero,
# such variants were not put into the 'smutations'
# field of a genome, due to the way the fwdpp mutation
# base class constructor was being called.  In PR #433,
# we now require a DFE to tell the Mutation type if the
# variant is neutral or not, and the expectation is that
# the answer is "False".  This test creates a new Sregion
# type where all variants have effect size zero and we run
# a quick sim to make sure that they are all in both genomes
# and in the mutation table.

import unittest

import EsizeZero
import fwdpy11


class Recorder(object):
    def __init__(self):
        self.data = []

    def __call__(self, pop, sampler):
        in_genome = [False] * len(pop.mutations)
        for i in pop.diploids:
            for g in [i.first, i.second]:
                for k in pop.haploid_genomes[g].smutations:
                    in_genome[k] = True
        temp = [in_genome[i.key] for i in pop.tables.mutations]
        self.data.append(all([i is True for i in temp]) is True)


class TestEsizeZero(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        pdict = {
            "nregions": [],
            "sregions": [EsizeZero.EsizeZero(0, 1, 1, True, 0)],
            "recregions": [fwdpy11.PoissonInterval(0, 1, 1)],
            "rates": (0, 5, None),
            "gvalue": fwdpy11.Multiplicative(2.0),
            "demography": fwdpy11.DiscreteDemography(),
            "simlen": 10,
        }
        self.params = fwdpy11.ModelParams(**pdict)
        self.pop = fwdpy11.DiploidPopulation(500, 1.0)
        self.rng = fwdpy11.GSLrng(512035)

    def test_simulation(self):
        r = Recorder()
        fwdpy11.evolvets(self.rng, self.pop, self.params, 1, r)
        self.assertTrue(all([i is True for i in r.data]))


if __name__ == "__main__":
    unittest.main()
