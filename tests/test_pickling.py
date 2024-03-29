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


class TestPickleDiploidPopulation(unittest.TestCase):
    @classmethod
    def setUp(self):
        from quick_pops import quick_neutral_slocus

        self.pop = quick_neutral_slocus()

    def testPickleDiploidsPy(self):
        import pickle

        pd = pickle.dumps(self.pop.diploids)
        p = pickle.loads(pd)
        for i, j in zip(p, self.pop.diploids):
            self.assertEqual(i, j)

    def testPicklePopPy(self):
        import pickle

        p = pickle.dumps(self.pop, -1)
        pp = pickle.loads(p)
        self.assertEqual(pp, self.pop)

    def testPickleTableCollection(self):
        import pickle

        p = pickle.dumps(self.pop.tables)
        up = pickle.loads(p)
        self.assertEqual(up.genome_length, self.pop.tables.genome_length)


class TestPickleDiploidPopulationTreeSequences(unittest.TestCase):
    @classmethod
    def setUp(self):
        import fwdpy11

        self.N = 1000
        self.rho = 1.0
        self.theta = 100.0
        self.nreps = 500
        self.mu = self.theta / (4 * self.N)
        self.r = self.rho / (4 * self.N)

        self.GSS = fwdpy11.GaussianStabilizingSelection.single_trait(
            [fwdpy11.Optimum(VS=1.0, optimum=0.0, when=0)]
        )
        a = fwdpy11.Additive(2.0, self.GSS)
        demography = fwdpy11.ForwardDemesGraph.tubes(
            [self.N], burnin=100, burnin_is_exact=True
        )
        self.p = {
            "nregions": [],
            "sregions": [fwdpy11.GaussianS(0, 1, 1, 0.25)],
            "recregions": [fwdpy11.PoissonInterval(0, 1, self.r)],
            "rates": (0.0, 0.025, None),
            "gvalue": a,
            "prune_selected": False,
            "demography": demography,
            "simlen": demography.final_generation,
        }
        self.params = fwdpy11.ModelParams(**self.p)
        self.rng = fwdpy11.GSLrng(101 * 45 * 110 * 210)
        self.pop = fwdpy11.DiploidPopulation(self.N, 1.0)
        fwdpy11.evolvets(self.rng, self.pop, self.params, 100)

    def testPickleNodeTable(self):
        import pickle
        import numpy as np

        pn = pickle.dumps(self.pop.tables.nodes)
        up = pickle.loads(pn)
        nodes = np.array(self.pop.tables.nodes)
        nodes2 = np.array(up)
        self.assertTrue(np.array_equal(nodes, nodes2))

    def testPickleEdgeTable(self):
        import pickle
        import numpy as np

        p = pickle.dumps(self.pop.tables.edges)
        u = pickle.loads(p)
        edges = np.array(self.pop.tables.edges)
        edges2 = np.array(u)
        self.assertTrue(np.array_equal(edges, edges2))

    def testPickleMutationTable(self):
        import pickle
        import numpy as np

        p = pickle.dumps(self.pop.tables.mutations)
        u = pickle.loads(p)
        mutations = np.array(self.pop.tables.mutations)
        mutations2 = np.array(u)
        self.assertTrue(np.array_equal(mutations, mutations2))

    def testPickleTableCollection(self):
        import pickle

        p = pickle.dumps(self.pop.tables)
        u = pickle.loads(p)
        self.assertTrue(u == self.pop.tables)

    def testPickleMetaData(self):
        import pickle
        import numpy as np

        p = pickle.dumps(self.pop.diploid_metadata)
        up = pickle.loads(p)
        md = np.array(self.pop.diploid_metadata)
        md2 = np.array(up)
        self.assertTrue(np.array_equal(md, md2))

    def testPicklingToFile(self):
        import fwdpy11
        import lzma
        import os

        fname = "pickleFileTesting.pickle"
        with lzma.open(fname, "wb") as f:
            self.pop.pickle_to_file(f)
        with lzma.open(fname, "rb") as f:
            pop = fwdpy11.DiploidPopulation.load_from_pickle_file(f)
        self.assertEqual(pop.N, self.pop.N)
        self.assertEqual(pop.generation, self.pop.generation)
        self.assertTrue(pop.diploids == self.pop.diploids)
        self.assertTrue(pop.haploid_genomes == self.pop.haploid_genomes)
        self.assertTrue(pop.mutations == self.pop.mutations)
        self.assertTrue(np.array_equal(pop.mcounts, self.pop.mcounts))
        self.assertTrue(pop.tables == self.pop.tables)
        self.assertTrue(pop == self.pop)
        os.remove(fname)


if __name__ == "__main__":
    unittest.main()
