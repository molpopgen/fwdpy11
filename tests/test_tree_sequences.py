import unittest
import fwdpy11
import fwdpy11.genetic_values
import fwdpy11.ts
import fwdpy11.model_params
import fwdpy11.wright_fisher_ts
import msprime
import numpy as np


class testTreeSequences(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.N = 1000
        self.demography = np.array([self.N]*25, dtype=np.uint32)
        self.rho = 100.
        self.theta = 100.
        self.nreps = 500
        self.mu = self.theta/(4*self.N)
        self.r = self.rho/(4*self.N)

        self.GSS = fwdpy11.genetic_values.GSS(VS=1, opt=0)
        a = fwdpy11.genetic_values.SlocusAdditive(2.0, self.GSS)
        self.p = {'nregions': [],
                  'sregions': [fwdpy11.GaussianS(0, 1, 1, 0.25)],
                  'recregions': [fwdpy11.Region(0, 1, 1)],
                  'rates': (0.0, 0.0, self.r),
                  'gvalue': a,
                  'prune_selected': False,
                  'demography': self.demography
                  }
        self.params = fwdpy11.model_params.ModelParams(**self.p)
        self.rng = fwdpy11.GSLrng(101*45*110*210)
        self.pop = fwdpy11.SlocusPop(self.N, 1.0)
        fwdpy11.wright_fisher_ts.evolve(self.rng, self.pop, self.params, 100)

    def test_dump_to_msprime(self):
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
        nview = np.array(self.pop.tables.nodes, copy=False)
        st1 = np.sort(nview['time'])
        st2 = np.sort(dumped_ts.tables.nodes.time)
        self.assertTrue(np.array_equal(st1, st2))

    def test_simplify_to_sample(self):
        dumped_ts = self.pop.dump_tables_to_msprime()
        tt = 0.0
        for i in self.pop.tables.nodes:
            tt += i.time
        print(tt)
        print(dumped_ts.tables.nodes.time.sum())
        tc = msprime.TableCollection(self.pop.tables.genome_length())
        tc.nodes.set_columns(flags=dumped_ts.tables.nodes.flags,
                             time=dumped_ts.tables.nodes.time)
        tc.edges.set_columns(parent=dumped_ts.tables.edges.parent,
                             child=dumped_ts.tables.edges.child,
                             left=dumped_ts.tables.edges.left,
                             right=dumped_ts.tables.edges.right)
        tc.sort()
        ts = tc.tree_sequence()
        nv = np.array(self.pop.tables.nodes, copy=False)
        print(ts.tables.nodes.time.sum(), nv['time'].sum())
        samples = np.arange(0, 2*self.pop.N, 50, dtype=np.int32)
        mspts = ts.simplify(samples=samples.tolist())
        fp11ts, maps = fwdpy11.ts.simplify(self.pop, samples)
        nv = np.array(fp11ts.nodes, copy=False)
        print(nv)
        print("??", len(mspts.tables.nodes.time))
        print(mspts.tables.nodes.time.sum(), nv['time'].sum())


if __name__ == "__main__":
    unittest.main()
