import unittest
import fwdpy11
import fwdpy11.genetic_values
import fwdpy11.ts
import fwdpy11.model_params
import fwdpy11.wright_fisher_ts
import numpy as np


class testTreeSequences(unittest.TestCase):
    @classmethod
    def setUp(self):
        # TODO add neutral variants
        self.N = 1000
        self.demography = np.array([self.N]*self.N, dtype=np.uint32)
        self.rho = 1.
        self.theta = 100.
        self.nreps = 500
        self.mu = self.theta/(4*self.N)
        self.r = self.rho/(4*self.N)

        self.GSS = fwdpy11.genetic_values.GSS(VS=1, opt=0)
        a = fwdpy11.genetic_values.SlocusAdditive(2.0, self.GSS)
        self.p = {'nregions': [],
                  'sregions': [fwdpy11.GaussianS(0, 1, 1, 0.25)],
                  'recregions': [fwdpy11.Region(0, 1, 1)],
                  'rates': (0.0, 0.025, self.r),
                  'gvalue': a,
                  'prune_selected': False,
                  'demography': self.demography
                  }
        self.params = fwdpy11.model_params.ModelParams(**self.p)
        self.rng = fwdpy11.GSLrng(101*45*110*210)
        self.pop = fwdpy11.SlocusPop(self.N, 1.0)
        fwdpy11.wright_fisher_ts.evolve(self.rng, self.pop, self.params, 100)

    def test_dump_to_msprime(self):
        # TODO: test leaf counts of mutations in msprmie
        # vs fwdpy11 and cross-references with self.pop.mcounts
        dumped_ts = self.pop.dump_tables_to_msprime()
        self.assertEqual(len(dumped_ts.tables.nodes),
                         len(self.pop.tables.nodes))
        self.assertEqual(len(dumped_ts.tables.edges),
                         len(self.pop.tables.edges))
        self.assertEqual(len(dumped_ts.tables.mutations),
                         len(self.pop.tables.mutations))
        eview = np.array(self.pop.tables.edges, copy=False)
        self.assertEqual(eview['parent'].sum(),
                         dumped_ts.tables.edges.parent.sum())
        self.assertEqual(eview['child'].sum(),
                         dumped_ts.tables.edges.child.sum())
        self.assertEqual(eview['left'].sum(),
                         dumped_ts.tables.edges.left.sum())
        self.assertEqual(eview['right'].sum(),
                         dumped_ts.tables.edges.right.sum())
        tv = fwdpy11.ts.TreeVisitor(self.pop.tables,
                                    [i for i in range(2*self.pop.N)])
        tt_fwd = 0
        while tv(False) is True:
            m = tv.tree()
            tt_fwd += m.total_time(self.pop.tables.nodes)
        tt_msprime = 0
        for t in dumped_ts.trees():
            tt_msprime += t.get_total_branch_length()
        self.assertEqual(tt_fwd, tt_msprime)

    def test_leaf_counts_vs_mcounts(self):
        tv = fwdpy11.ts.TreeVisitor(self.pop.tables,
                                    [i for i in range(2*self.pop.N)])
        mv = np.array(self.pop.tables.mutations, copy=False)
        muts = np.array(self.pop.mutations.array())
        p = muts['pos']
        while tv(False) is True:
            m = tv.tree()
            l, r = m.left, m.right
            mt = [i for i in mv if p[i[1]] >= l and p[i[1]] < r]
            for i in mt:
                self.assertEqual(m.leaf_counts[i[0]],
                                 self.pop.mcounts[i[1]])

    def test_simplify_to_sample(self):
        """
        Simplify to a sample using fwdpy11.ts and
        msprime, then test that total time on output
        is the same from both sources and that
        the mutation tables contain the same
        positions after simplification.
        """
        dumped_ts = self.pop.dump_tables_to_msprime()
        tt = 0.0
        for i in self.pop.tables.nodes:
            tt += i.time
        samples = np.arange(0, 2*self.pop.N, 50, dtype=np.int32)
        mspts = dumped_ts.simplify(samples=samples.tolist())
        fp11ts, idmap = fwdpy11.ts.simplify(self.pop, samples)
        for i in range(len(fp11ts.edges)):
            self.assertTrue(fp11ts.edges[i].parent < len(fp11ts.nodes))
            self.assertTrue(fp11ts.edges[i].child < len(fp11ts.nodes))
        for s in samples:
            self.assertEqual(
                fp11ts.nodes[idmap[s]].time, self.pop.generation)
        tt_fwd = 0.0
        tv = fwdpy11.ts.TreeVisitor(fp11ts, [i for i in range(len(samples))])
        while tv(False) is True:
            m = tv.tree()
            tt_fwd += m.total_time(fp11ts.nodes)
        tt_msprime = 0.0
        for t in mspts.trees():
            tt_msprime += t.get_total_branch_length()
        self.assertEqual(tt_fwd, tt_msprime)

        self.assertEqual(len(fp11ts.mutations),
                         len(mspts.tables.mutations))
        fp11_pos = np.array([self.pop.mutations[i.key].pos
                             for i in fp11ts.mutations])
        fp11_pos = np.sort(fp11_pos)
        msp_pos = np.sort(mspts.tables.sites.position)
        self.assertTrue(np.array_equal(fp11_pos, msp_pos))

    def test_genotype_matrix(self):
        """
        Make data matrix objects from the tree sequences
        and compare their contents to those of msprime
        as well as to an explicit calculation of mutation counts.
        """
        dm = fwdpy11.ts.make_data_matrix(self.pop,
                                         [i for i in range(2*self.pop.N)],
                                         False, True)
        sa = np.array(dm.selected)
        cs = np.sum(sa, axis=1)
        dumped_ts = self.pop.dump_tables_to_msprime()
        mm = dumped_ts.genotype_matrix()
        mc = np.sum(mm, axis=1)
        ec = np.zeros(len(self.pop.mutations), dtype=np.uint32)
        for d in self.pop.diploids:
            for m in self.pop.gametes[d.first].smutations:
                ec[m] += 1
            for m in self.pop.gametes[d.second].smutations:
                ec[m] += 1
        for i, j, k in zip(self.pop.tables.mutations, cs, mc):
            self.assertEqual(self.pop.mcounts[i.key],
                             ec[i.key])
            self.assertEqual(ec[i.key], j)
            self.assertEqual(ec[i.key], k)
        self.assertTrue(np.array_equal(ec, np.array(self.pop.mcounts)))
        self.assertEqual(ec.sum(), cs.sum())
        self.assertEqual(ec.sum(), mc.sum())
        self.assertTrue(np.array_equal(sa, mm))
        self.assertTrue(np.array_equal(cs, mc))

    def test_VariantIterator(self):
        """
        Test VariantIterator by asserting
        that sum of genotypes equal values in
        the corresponding DataMatrix and
        those in pop.mcounts
        """
        dm = fwdpy11.ts.make_data_matrix(self.pop,
                                         [i for i in range(2*self.pop.N)],
                                         False, True)
        sa = np.array(dm.selected)
        cs = np.sum(sa, axis=1)
        i = 0
        vi = fwdpy11.ts.VariantIterator(self.pop.tables,
                                        self.pop.mutations,
                                        [i for i in range(2*self.pop.N)])
        for v in vi:
            c = self.pop.mcounts[self.pop.tables.mutations[i].key]
            self.assertEqual(c, cs[i])
            self.assertEqual(c, v.genotypes.sum())
            i += 1
        self.assertEqual(i, len(self.pop.tables.mutations))

    def test_VariantIteratorFromPopulation(self):
        dm = fwdpy11.ts.make_data_matrix(self.pop,
                                         [i for i in range(2*self.pop.N)],
                                         False, True)
        sa = np.array(dm.selected)
        cs = np.sum(sa, axis=1)
        i = 0
        vi = fwdpy11.ts.VariantIterator(self.pop)
        for v in vi:
            c = self.pop.mcounts[self.pop.tables.mutations[i].key]
            self.assertEqual(c, cs[i])
            self.assertEqual(c, v.genotypes.sum())
            i += 1
        self.assertEqual(i, len(self.pop.tables.mutations))


if __name__ == "__main__":
    unittest.main()
