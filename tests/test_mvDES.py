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


def gvalue_multiplicative(pop, ind, scaling):
    g = 1.0
    keys = [k for k in pop.haploid_genomes[pop.diploids[ind].first].smutations]
    keys.extend([k for k in pop.haploid_genomes[pop.diploids[ind].second].smutations])
    keys = np.array(keys, dtype=np.uint32)
    key_counts = np.unique(keys, return_counts=True)

    ind_md = pop.diploid_metadata[ind]
    for i, j in zip(key_counts[0], key_counts[1]):
        if j == 1:
            g *= (1.0 + pop.mutations[i].esizes[ind_md.deme]
                  * pop.mutations[i].heffects[ind_md.deme])
        elif j == 2:
            g *= (1.0 + scaling*pop.mutations[i].esizes[ind_md.deme])

    return g


class TestMultiplicativeWithExpSNoMigration(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.pop = fwdpy11.DiploidPopulation([50, 50], 1)
        mvDES = fwdpy11.mvDES([fwdpy11.ExpS(0, 1, 1, 0.1), fwdpy11.ExpS(
            0, 1, 1, -0.1)], np.zeros(2), np.identity(2))

        self.pdict = {'nregions': [],
                      'sregions': [mvDES],
                      'recregions': [fwdpy11.PoissonInterval(0, 1, 1e-2)],
                      'rates': (0, 1e-2, None),
                      'gvalue': fwdpy11.Multiplicative(2., ndemes=2),
                      'demography': fwdpy11.DiscreteDemography(),
                      'simlen': 100,
                      'prune_selected': True}
        self.params = fwdpy11.ModelParams(**self.pdict)
        self.rng = fwdpy11.GSLrng(918273)
        fwdpy11.evolvets(self.rng, self.pop, self.params, 100)
        n0 = 0
        n1 = 0
        for i, md in enumerate(self.pop.diploid_metadata):
            ng0 = len(self.pop.haploid_genomes[self.pop.diploids[i].first].smutations)
            ng1 = len(self.pop.haploid_genomes[self.pop.diploids[i].second].smutations)
            if ng0 + ng1 > 0:
                if md.deme == 0:
                    n0 += 1
                elif md.deme == 1:
                    n1 += 1
        if n0 == 0 or n1 == 0:
            self.assertFail("we don't have individuals with mutations in each deme")

    def test_it(self):
        gv = np.zeros(self.pop.N)
        for i in range(self.pop.N):
            gv[i] = gvalue_multiplicative(self.pop, i, 2.0)
        md = np.array(self.pop.diploid_metadata, copy=False)
        self.assertTrue(np.allclose(gv, md['g']))


if __name__ == "__main__":
    unittest.main()
