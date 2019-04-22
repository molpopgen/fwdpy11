import unittest
import fwdpy11
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

        self.GSS = fwdpy11.GSS(VS=1, opt=0)
        a = fwdpy11.Additive(2.0, self.GSS)
        self.p = {'nregions': [],
                  'sregions': [fwdpy11.GaussianS(0, 1, 1, 0.25)],
                  'recregions': [fwdpy11.Region(0, 1, 1)],
                  'rates': (0.0, 0.025, self.r),
                  'gvalue': a,
                  'prune_selected': False,
                  'demography': self.demography
                  }
        self.params = fwdpy11.ModelParams(**self.p)
        self.rng = fwdpy11.GSLrng(101*45*110*210)
        self.pop = fwdpy11.DiploidPopulation(self.N, 1.0)
        fwdpy11.evolvets(self.rng, self.pop, self.params, 100)

    def test_simplify(self):
        tables, idmap = fwdpy11.simplify(self.pop, [i for i in range(10)])
        for i in range(10):
            self.assertTrue(idmap[i] != fwdpy11.NULL_NODE)

    def test_simplify_tables(self):
        tables, idmap = fwdpy11.simplify_tables(
            self.pop.tables, self.pop.mutations, [i for i in range(10)])
        for i in range(10):
            self.assertTrue(idmap[i] != fwdpy11.NULL_NODE)

    def test_simplify_numpy_array(self):
        tables, idmap = fwdpy11.simplify(
            self.pop, np.array([i for i in range(10)]))
        for i in range(10):
            self.assertTrue(idmap[i] != fwdpy11.NULL_NODE)

    def test_simplify_tables_numpy_array(self):
        tables, idmap = fwdpy11.simplify_tables(
            self.pop.tables, self.pop.mutations,
            np.array([i for i in range(10)]))
        for i in range(10):
            self.assertTrue(idmap[i] != fwdpy11.NULL_NODE)

    def test_dump_to_tskit(self):
        # TODO: test leaf counts of mutations in msprmie
        # vs fwdpy11 and cross-references with self.pop.mcounts
        dumped_ts = self.pop.dump_tables_to_tskit()
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
        tv = fwdpy11.TreeIterator(self.pop.tables,
                                  [i for i in range(2*self.pop.N)])
        tt_fwd = 0
        for t in tv:
            tt_fwd += t.total_time(self.pop.tables.nodes)
        tt_tskit = 0
        for t in dumped_ts.trees():
            tt_tskit += t.get_total_branch_length()
        self.assertEqual(tt_fwd, tt_tskit)

    def test_TreeIterator(self):
        with self.assertRaises(ValueError):
            tv = fwdpy11.TreeIterator(self.pop.tables,
                                      [i for i in range(2*self.pop.N)],
                                      True, 1, 0)
        with self.assertRaises(ValueError):
            tv = fwdpy11.TreeIterator(self.pop.tables,
                                      [i for i in range(2*self.pop.N)],
                                      False, 1, 0)

        for i in np.arange(0., 1., 0.1):
            tv = fwdpy11.TreeIterator(self.pop.tables,
                                      [i for i in range(2*self.pop.N)],
                                      True, i, i+0.1)
            for ti in tv:
                a = ti.left < i+0.1
                b = i < ti.right
                self.assertTrue(a and b)

            tv = fwdpy11.TreeIterator(self.pop.tables,
                                      [i for i in range(2*self.pop.N)],
                                      False, i, i+0.1)
            for ti in tv:
                a = ti.left < i+0.1
                b = i < ti.right
                self.assertTrue(a and b)

    def test_leaf_counts_vs_mcounts(self):
        tv = fwdpy11.TreeIterator(self.pop.tables,
                                  [i for i in range(2*self.pop.N)])
        mv = np.array(self.pop.tables.mutations, copy=False)
        muts = self.pop.mutations_ndarray
        p = muts['pos']
        for t in tv:
            l, r = t.left, t.right
            mt = [i for i in mv if p[i[1]] >= l and p[i[1]] < r]
            for i in mt:
                self.assertEqual(t.leaf_counts(i[0]),
                                 self.pop.mcounts[i[1]])

    def test_simplify_to_sample(self):
        """
        Simplify to a sample using fwdpy11 and
        tskit, then test that total time on output
        is the same from both sources and that
        the mutation tables contain the same
        positions after simplification.
        """
        dumped_ts = self.pop.dump_tables_to_tskit()
        tt = 0.0
        for i in self.pop.tables.nodes:
            tt += i.time
        samples = np.arange(0, 2*self.pop.N, 50, dtype=np.int32)
        mspts = dumped_ts.simplify(samples=samples.tolist())
        fp11ts, idmap = fwdpy11.simplify(self.pop, samples)
        for i in range(len(fp11ts.edges)):
            self.assertTrue(fp11ts.edges[i].parent < len(fp11ts.nodes))
            self.assertTrue(fp11ts.edges[i].child < len(fp11ts.nodes))
        for s in samples:
            self.assertEqual(
                fp11ts.nodes[idmap[s]].time, self.pop.generation)
        tt_fwd = 0.0
        tv = fwdpy11.TreeIterator(fp11ts, [i for i in range(len(samples))])
        for t in tv:
            tt_fwd += t.total_time(fp11ts.nodes)
        tt_tskit = 0.0
        for t in mspts.trees():
            tt_tskit += t.get_total_branch_length()
        self.assertEqual(tt_fwd, tt_tskit)

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
        and compare their contents to those of tskit
        as well as to an explicit calculation of mutation counts.
        """
        dm = fwdpy11.make_data_matrix(self.pop,
                                      [i for i in range(2*self.pop.N)],
                                      False, True)
        sa = np.array(dm.selected)
        cs = np.sum(sa, axis=1)
        dumped_ts = self.pop.dump_tables_to_tskit()
        mm = dumped_ts.genotype_matrix()
        mc = np.sum(mm, axis=1)
        ec = np.zeros(len(self.pop.mutations), dtype=np.uint32)
        for d in self.pop.diploids:
            for m in self.pop.haploid_genomes[d.first].smutations:
                ec[m] += 1
            for m in self.pop.haploid_genomes[d.second].smutations:
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

    def test_genotype_matrix_ranges(self):
        dm = fwdpy11.data_matrix_from_tables(self.pop.tables,
                                             self.pop.mutations,
                                             [i for i in range(
                                                 2*self.pop.N)],
                                             False, True)
        spos = np.array(dm.selected.positions)
        for i in np.arange(0, self.pop.tables.genome_length, 0.1):
            dmi = fwdpy11.data_matrix_from_tables(self.pop.tables,
                                                  self.pop.mutations,
                                                  [i for i in range(
                                                      2*self.pop.N)],
                                                  record_neutral=False,
                                                  record_selected=True,
                                                  begin=i, end=i+0.1)
            w = np.where((spos >= i) & (spos < i+0.1))[0]
            self.assertTrue(np.array_equal(
                spos[w], np.array(dmi.selected.positions)))

    def test_VariantIterator(self):
        """
        Test VariantIterator by asserting
        that sum of genotypes equal values in
        the corresponding DataMatrix and
        those in pop.mcounts
        """
        dm = fwdpy11.data_matrix_from_tables(self.pop.tables,
                                             self.pop.mutations,
                                             [i for i in range(
                                                 2*self.pop.N)],
                                             False, True)
        sa = np.array(dm.selected)
        cs = np.sum(sa, axis=1)
        i = 0
        vi = fwdpy11.VariantIterator(self.pop.tables,
                                     self.pop.mutations,
                                     [i for i in range(2*self.pop.N)])
        for v in vi:
            c = self.pop.mcounts[self.pop.tables.mutations[i].key]
            self.assertEqual(c, cs[i])
            self.assertEqual(c, v.genotypes.sum())
            i += 1
        mc = np.array(self.pop.mcounts)
        self.assertEqual(i, len(np.where(mc > 0)[0]))
        self.assertEqual(i, len(self.pop.tables.mutations))

    def test_VariantIteratorFromPopulation(self):
        dm = fwdpy11.data_matrix_from_tables(self.pop.tables,
                                             self.pop.mutations,
                                             [i for i in range(
                                                 2*self.pop.N)],
                                             False, True)
        sa = np.array(dm.selected)
        cs = np.sum(sa, axis=1)
        i = 0
        vi = fwdpy11.VariantIterator(self.pop)
        for v in vi:
            c = self.pop.mcounts[self.pop.tables.mutations[i].key]
            self.assertEqual(c, cs[i])
            self.assertEqual(c, v.genotypes.sum())
            i += 1
        mc = np.array(self.pop.mcounts)
        self.assertEqual(i, len(np.where(mc > 0)[0]))
        self.assertEqual(i, len(self.pop.tables.mutations))

    def test_VariantIteratorBeginEnd(self):
        for i in np.arange(0, self.pop.tables.genome_length, 0.1):
            vi = fwdpy11.VariantIterator(self.pop.tables, self.pop.mutations,
                                         [i for i in range(2*self.pop.N)], i, i+0.1)
            nm = len([j for j in self.pop.tables.mutations if self.pop.mutations[j.key].pos >= i and
                      self.pop.mutations[j.key].pos < i+0.1])
            nseen = 0
            for v in vi:
                r = v.record
                self.assertTrue(self.pop.mutations[r.key].pos >= i)
                self.assertTrue(self.pop.mutations[r.key].pos < i+0.1)
                nseen += 1
            self.assertEqual(nm, nseen)

        # test bad start/stop
        with self.assertRaises(ValueError):
            vi = fwdpy11.VariantIterator(self.pop.tables, self.pop.mutations,
                                         [i for i in range(2*self.pop.N)], begin=0.5, end=0.25)

    def test_count_mutations(self):
        mc = fwdpy11.count_mutations(self.pop,
                                     [i for i in range(2*self.pop.N)])
        pmc = np.array(self.pop.mcounts)
        self.assertTrue(np.array_equal(mc, pmc))


class testSamplePreservation(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.N = 1000
        self.demography = np.array([self.N]*100, dtype=np.uint32)
        self.rho = 1.
        self.theta = 100.
        self.nreps = 500
        self.mu = self.theta/(4*self.N)
        self.r = self.rho/(4*self.N)
        self.GSS = fwdpy11.GSS(VS=1, opt=0)
        a = fwdpy11.Additive(2.0, self.GSS)
        self.p = {'nregions': [],
                  'sregions': [fwdpy11.GaussianS(0, 1, 1, 0.25)],
                  'recregions': [fwdpy11.Region(0, 1, 1)],
                  'rates': (0.0, 0.025, self.r),
                  'gvalue': a,
                  'prune_selected': False,
                  'demography': self.demography
                  }
        self.params = fwdpy11.ModelParams(**self.p)
        self.rng = fwdpy11.GSLrng(101*45*110*210)
        self.pop = fwdpy11.DiploidPopulation(self.N, 1.0)
        self.recorder = fwdpy11.RandomAncientSamples(seed=42,
                                                     samplesize=10,
                                                     timepoints=[i for i in range(1, 101)])
        fwdpy11.evolvets(
            self.rng, self.pop, self.params, 100, self.recorder)

    def test_Simulation(self):
        self.assertEqual(self.pop.generation, 100)

    def test_count_mutations_preserved_samples(self):
        mc = fwdpy11.count_mutations(self.pop,
                                     self.pop.tables.preserved_nodes)
        pmc = np.array(self.pop.mcounts_ancient_samples)
        self.assertTrue(np.array_equal(mc, pmc))

    def test_VariantIteratorFromPreservedSamples(self):
        n = np.array(self.pop.tables.nodes)
        pn = np.array(self.pop.tables.preserved_nodes)
        at = n['time'][pn]
        for u in np.unique(at):
            n = pn[np.where(at == u)[0]]
            vi = fwdpy11.VariantIterator(self.pop.tables,
                                         self.pop.mutations, n)
            for variant in vi:
                k = variant.record
                self.assertNotEqual(k.node, fwdpy11.NULL_NODE)
                self.assertNotEqual(k.key, np.iinfo(np.uint64).max)


class testSimplificationInterval(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.N = 1000
        self.demography = np.array([self.N]*100, dtype=np.uint32)
        self.rho = 1.
        self.theta = 100.
        self.nreps = 500
        self.mu = self.theta/(4*self.N)
        self.r = self.rho/(4*self.N)
        self.GSS = fwdpy11.GSS(VS=1, opt=0)
        a = fwdpy11.Additive(2.0, self.GSS)
        self.p = {'nregions': [],
                  'sregions': [fwdpy11.GaussianS(0, 1, 1, 0.25)],
                  'recregions': [fwdpy11.Region(0, 1, 1)],
                  'rates': (0.0, 0.025, self.r),
                  'gvalue': a,
                  'prune_selected': False,
                  'demography': self.demography
                  }
        self.params = fwdpy11.ModelParams(**self.p)
        self.rng = fwdpy11.GSLrng(101*45*110*210)
        self.pop = fwdpy11.DiploidPopulation(self.N, 1.0)
        self.recorder = fwdpy11.RandomAncientSamples(seed=42,
                                                     samplesize=10,
                                                     timepoints=[i for i in range(1, 101)])

    def testEvolve(self):
        # TODO: actually test something here :)
        fwdpy11.evolvets(
            self.rng, self.pop, self.params, 1, self.recorder)
        samples = [i for i in range(2*self.pop.N)] + \
            self.pop.tables.preserved_nodes
        vi = fwdpy11.TreeIterator(self.pop.tables, samples)


class testFixationPreservation(unittest.TestCase):
    def testQtraitSim(self):
        N = 1000
        demography = np.array([N]*10*N, dtype=np.uint32)
        rho = 1.
        r = rho/(4*N)

        GSS = fwdpy11.GSS(VS=1, opt=1)
        a = fwdpy11.Additive(2.0, GSS)
        p = {'nregions': [],
             'sregions': [fwdpy11.GaussianS(0, 1, 1, 0.25)],
             'recregions': [fwdpy11.Region(0, 1, 1)],
             'rates': (0.0, 0.005, r),
             'gvalue': a,
             'prune_selected': False,
             'demography': demography
             }
        params = fwdpy11.ModelParams(**p)
        rng = fwdpy11.GSLrng(101*45*110*210)
        pop = fwdpy11.DiploidPopulation(N, 1.0)
        fwdpy11.evolvets(rng, pop, params, 100,
                         track_mutation_counts=True)
        mc = fwdpy11.count_mutations(pop.tables, pop.mutations,
                                     [i for i in range(2*pop.N)])
        assert len(pop.fixations) > 0, "Test is meaningless without fixations"
        fixations = np.where(mc == 2*pop.N)[0]
        self.assertEqual(len(fixations), len(pop.fixations))

        # Brute-force calculation of fixations
        brute_force = np.zeros(len(pop.mutations), dtype=np.int32)
        for g in pop.haploid_genomes:
            if g.n > 0:
                for k in g.smutations:
                    brute_force[k] += g.n

        self.assertTrue(np.array_equal(brute_force, mc))
        self.assertTrue(np.array_equal(brute_force, pop.mcounts))

    def testPopGenSim(self):
        N = 1000
        demography = np.array([N]*10*N, dtype=np.uint32)
        rho = 1.
        r = rho/(4*N)

        a = fwdpy11.Multiplicative(2.0)
        p = {'nregions': [],
             'sregions': [fwdpy11.ExpS(0, 1, 1, 0.01)],
             'recregions': [fwdpy11.Region(0, 1, 1)],
             'rates': (0.0, 0.00005, r),
             'gvalue': a,
             'prune_selected': True,
             'demography': demography
             }
        params = fwdpy11.ModelParams(**p)
        rng = fwdpy11.GSLrng(101*45*110*210)
        pop = fwdpy11.DiploidPopulation(N, 1.0)
        fwdpy11.evolvets(rng, pop, params, 100,
                         track_mutation_counts=True)
        mc = fwdpy11.count_mutations(pop.tables, pop.mutations,
                                     [i for i in range(2*pop.N)])
        assert len(pop.fixations) > 0, "Test is meaningless without fixations"
        fixations = np.where(mc == 2*pop.N)[0]
        self.assertEqual(len(fixations), 0)

        # Brute-force calculation of fixations
        brute_force = np.zeros(len(pop.mutations), dtype=np.int32)
        for g in pop.haploid_genomes:
            if g.n > 0:
                for k in g.smutations:
                    brute_force[k] += g.n

        self.assertTrue(np.array_equal(brute_force, mc))
        self.assertTrue(np.array_equal(brute_force, pop.mcounts))


class testMetaData(unittest.TestCase):
    """
    Integration test.

    Rebuild genetic values from mutations
    stored in tree sequence
    """

    def testQtraitSim(self):
        N = 1000
        demography = np.array([N]*10*N, dtype=np.uint32)
        rho = 1.
        r = rho/(4*N)

        GSS = fwdpy11.GSS(VS=1, opt=1)
        a = fwdpy11.Additive(2.0, GSS)
        p = {'nregions': [],
             'sregions': [fwdpy11.GaussianS(0, 1, 1, 0.25)],
             'recregions': [fwdpy11.PoissonInterval(0, 1, r)],
             'rates': (0.0, 0.005, None),
             'gvalue': a,
             'prune_selected': False,
             'demography': demography
             }
        params = fwdpy11.ModelParams(**p)
        rng = fwdpy11.GSLrng(101*45*110*210)
        pop = fwdpy11.DiploidPopulation(N, 1.0)

        class Recorder(object):
            """ Records entire pop every 100 generations """

            def __call__(self, pop, recorder):
                if pop.generation % 100 == 0.0:
                    recorder.assign(np.arange(pop.N, dtype=np.int32))

        r = Recorder()
        fwdpy11.evolvets(rng, pop, params, 100, r)

        ancient_sample_metadata = np.array(
            pop.ancient_sample_metadata, copy=False)
        alive_sample_metadata = np.array(pop.diploid_metadata, copy=False)
        metadata = np.hstack((ancient_sample_metadata, alive_sample_metadata))

        nodes = np.array(pop.tables.nodes, copy=False)
        metadata_nodes = metadata['nodes'].flatten()
        metadata_node_times = nodes['time'][metadata_nodes]
        metadata_record_times = nodes['time'][metadata['nodes'][:, 0]]

        genetic_trait_values_from_sim = []
        genetic_values_from_ts = []
        for u in np.unique(metadata_node_times):
            samples_at_time_u = metadata_nodes[np.where(
                metadata_node_times == u)]
            vi = fwdpy11.VariantIterator(
                pop.tables, pop.mutations, samples_at_time_u)
            sum_esizes = np.zeros(len(samples_at_time_u))
            for variant in vi:
                g = variant.genotypes
                r = variant.record
                mutant = np.where(g == 1)[0]
                sum_esizes[mutant] += pop.mutations[r.key].s
            ind = int(len(samples_at_time_u)/2)
            temp_gvalues = np.zeros(ind)
            temp_gvalues += sum_esizes[0::2]
            temp_gvalues += sum_esizes[1::2]
            genetic_values_from_ts.extend(temp_gvalues.tolist())
            genetic_trait_values_from_sim.extend(
                metadata['g'][np.where(metadata_record_times == u)[0]].tolist())

        for i, j in zip(genetic_trait_values_from_sim, genetic_values_from_ts):
            self.assertAlmostEqual(i, j)


if __name__ == "__main__":
    unittest.main()
