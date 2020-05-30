#
# Copyright(C) 2020 Kevin Thornton < krthornt @uci.edu>
#
# This file is part of fwdpy11.
#
# fwdpy11 is free software : you can redistribute it and / or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# fwdpy11 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with fwdpy11.If not, see < http: //www.gnu.org/licenses/>.
#

"""
Test written while addressing the following GitHub issues:
388
389
390
"""

import pickle
import tempfile
import unittest

import numpy as np

import fwdpy11


def set_up_quant_trait_model():
    N = 1000
    rho = 1.0
    r = rho / (4 * N)
    Opt = fwdpy11.Optimum
    GSSmo = fwdpy11.GSSmo(
        [Opt(when=0, optimum=0.0, VS=1.0), Opt(when=10 * N, optimum=1.0, VS=1.0)]
    )
    a = fwdpy11.Additive(2.0, GSSmo)
    p = {
        "nregions": [],
        "sregions": [fwdpy11.GaussianS(0, 1, 1, 0.25)],
        "recregions": [fwdpy11.Region(0, 1, 1)],
        "rates": (0.0, 0.025, r),
        "gvalue": a,
        "prune_selected": False,
        "demography": fwdpy11.DiscreteDemography(),
        "simlen": 10 * N + 100,
    }
    params = fwdpy11.ModelParams(**p)
    rng = fwdpy11.GSLrng(101 * 45 * 110 * 210)
    pop = fwdpy11.DiploidPopulation(N, 1.0)
    return params, rng, pop


def set_up_two_trait_quant_trait_model():
    N = 1000
    rho = 1.0
    r = rho / (4 * N)

    optima = [
        fwdpy11.PleiotropicOptima(when=0, optima=np.zeros(2), VS=2.0),
        fwdpy11.PleiotropicOptima(
            when=10 * N, optima=np.array([np.sqrt(2.0), 0]), VS=2.0
        ),
    ]
    GSSmo = fwdpy11.MultivariateGSSmo(optima)
    a = fwdpy11.StrictAdditiveMultivariateEffects(2, 0, GSSmo)
    vcov = np.identity(2)
    np.fill_diagonal(vcov, 0.25)
    DES = fwdpy11.MultivariateGaussianEffects(0, 1, 1, vcov)
    p = {
        "nregions": [],
        "sregions": [DES],
        "recregions": [fwdpy11.Region(0, 1, 1)],
        "rates": (0.0, 0.025, r),
        "gvalue": a,
        "prune_selected": False,
        "demography": fwdpy11.DiscreteDemography(),
        "simlen": 10 * N + 100,
    }
    params = fwdpy11.ModelParams(**p)
    rng = fwdpy11.GSLrng(101 * 45 * 110 * 210)
    pop = fwdpy11.DiploidPopulation(N, 1.0)
    return params, rng, pop


class PreserveN(object):
    def __init__(self, start, n):
        self.start = start
        self.n = n

    def __call__(self, pop, recorder):
        if pop.generation >= self.start:
            recorder.assign(np.arange(self.n, dtype=np.int32))


class TestNoPleiotropy(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.params, self.rng, self.pop = set_up_quant_trait_model()
        preserver = PreserveN(10 * self.pop.N, 10)
        fwdpy11.evolvets(
            self.rng, self.pop, self.params, 100, preserver, record_gvalue_matrix=True
        )

    def test_alive_genetic_values(self):
        for i, j in zip(self.pop.genetic_values.flatten(), self.pop.diploid_metadata):
            self.assertEqual(i, j.g)

    def test_alive_genetic_value_reconstruction(self):
        ti = fwdpy11.TreeIterator(
            self.pop.tables, self.pop.alive_nodes, update_samples=True
        )
        gv = np.zeros(2 * self.pop.N)
        for t in ti:
            for m in t.mutations():
                for b in t.samples_below(m.node):
                    gv[b] += self.pop.mutations[m.key].s
        gv = gv.reshape((self.pop.N, 2))
        gv = np.sum(gv, axis=1)
        md = np.array(self.pop.diploid_metadata, copy=False)
        self.assertTrue(np.allclose(gv, md["g"]))

    def test_ancient_sample_genetic_values(self):
        for i, j in zip(
            self.pop.ancient_sample_genetic_values, self.pop.ancient_sample_metadata
        ):
            self.assertEqual(i, j.g)

    def test_ancient_sample_genetic_value_reconstruction(self):
        as_gv = self.pop.ancient_sample_genetic_values
        amd = np.array(self.pop.ancient_sample_metadata, copy=False)
        nt = np.array(self.pop.tables.nodes, copy=False)
        ancient_nodes = amd["nodes"][:, 0]
        ancient_node_times = nt["time"][ancient_nodes]
        for time, n, md in self.pop.sample_timepoints(include_alive=False):
            w = np.where(ancient_node_times == time)[0]
            gvslice = as_gv[w].flatten()
            ti = fwdpy11.TreeIterator(self.pop.tables, n, update_samples=True)
            gv = np.zeros(len(n))
            node_map = np.array([np.iinfo(np.int32).max] * len(nt), dtype=np.int32)
            for i, j in enumerate(n):
                node_map[j] = i
            for t in ti:
                for m in t.mutations():
                    for b in t.samples_below(m.node):
                        gv[node_map[b]] += self.pop.mutations[m.key].s
            gv = gv.reshape((len(w), 2))
            gv = np.sum(gv, axis=1)
            self.assertTrue(np.allclose(gv, gvslice))

    def test_pickling(self):
        pp = pickle.dumps(self.pop, -1)
        up = pickle.loads(pp)
        gv = self.pop.genetic_values
        upgv = up.genetic_values
        self.assertTrue(np.array_equal(gv, upgv))
        gv = self.pop.ancient_sample_genetic_values
        upgv = up.ancient_sample_genetic_values
        self.assertTrue(np.array_equal(gv, upgv))

    def test_pickle_to_file(self):
        with tempfile.NamedTemporaryFile(delete=False) as tfile:
            self.pop.pickle_to_file(tfile)
            tfile.seek(0)
            pop2 = fwdpy11.DiploidPopulation.load_from_pickle_file(tfile)
        gv = self.pop.genetic_values
        gv2 = pop2.genetic_values
        self.assertTrue(np.array_equal(gv, gv2))
        gv = self.pop.ancient_sample_genetic_values
        gv2 = pop2.ancient_sample_genetic_values
        self.assertTrue(np.array_equal(gv, gv2))


class TestTwoTraitsIsotropy(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        """
        Differs from previous test:
        we only preserve after the
        optimum shifts, which makes
        it simpler to verify fitnesses.
        """
        self.params, self.rng, self.pop = set_up_two_trait_quant_trait_model()
        preserver = PreserveN(10 * self.pop.N + 1, 10)
        fwdpy11.evolvets(
            self.rng, self.pop, self.params, 100, preserver, record_gvalue_matrix=True
        )
        self.zopt = np.sqrt(2.0)
        self.VS = 2.0

    def test_alive_genetic_values_focal_trait(self):
        gv = self.pop.genetic_values
        for i, j in zip(range(gv.shape[0]), self.pop.diploid_metadata):
            self.assertEqual(gv[i, 0], j.g)

    def test_alive_genetic_value_reconstruction(self):
        gv = self.pop.genetic_values
        ti = fwdpy11.TreeIterator(
            self.pop.tables, self.pop.alive_nodes, update_samples=True
        )
        gv0 = np.zeros(2 * self.pop.N)
        gv1 = np.zeros(2 * self.pop.N)
        for t in ti:
            for m in t.mutations():
                for b in t.samples_below(m.node):
                    gv0[b] += self.pop.mutations[m.key].esizes[0]
                    gv1[b] += self.pop.mutations[m.key].esizes[1]
        # Check trait 0
        gv0 = gv0.reshape((self.pop.N, 2))
        gv0 = np.sum(gv0, axis=1)
        self.assertTrue(np.allclose(gv[:, 0], gv0))
        # Check trait 1
        gv1 = gv1.reshape((self.pop.N, 2))
        gv1 = np.sum(gv1, axis=1)
        self.assertTrue(np.allclose(gv[:, 1], gv1))

        # Now, we check fitness
        md = np.array(self.pop.diploid_metadata, copy=False)
        for i, j, w in zip(gv0, gv1, md["w"]):
            d0 = np.power(i - self.zopt, 2.0)
            d1 = np.power(j - 0.0, 2.0)
            self.assertTrue(np.isclose(np.exp(-(d0 + d1) / (2.0 * self.VS)), w))

    def test_ancient_sample_genetic_values_focal_trait(self):
        gv = self.pop.ancient_sample_genetic_values
        for i, j in zip(range(gv.shape[0]), self.pop.ancient_sample_metadata):
            self.assertEqual(gv[i, 0], j.g)

    def test_ancient_sample_genetic_value_reconstruction(self):
        as_gv = self.pop.ancient_sample_genetic_values
        amd = np.array(self.pop.ancient_sample_metadata, copy=False)
        nt = np.array(self.pop.tables.nodes, copy=False)
        ancient_nodes = amd["nodes"][:, 0]
        ancient_node_times = nt["time"][ancient_nodes]
        for time, n, md in self.pop.sample_timepoints(include_alive=False):
            w = np.where(ancient_node_times == time)[0]
            gvslice = as_gv[w, :]
            ti = fwdpy11.TreeIterator(self.pop.tables, n, update_samples=True)
            gv0 = np.zeros(len(n))
            gv1 = np.zeros(len(n))
            node_map = np.array([np.iinfo(np.int32).max] * len(nt), dtype=np.int32)
            for i, j in enumerate(n):
                node_map[j] = i
            for t in ti:
                for m in t.mutations():
                    for b in t.samples_below(m.node):
                        gv0[node_map[b]] += self.pop.mutations[m.key].esizes[0]
                        gv1[node_map[b]] += self.pop.mutations[m.key].esizes[1]
            # Check trait 0
            gv0 = gv0.reshape((len(w), 2))
            gv0 = np.sum(gv0, axis=1)
            self.assertTrue(np.allclose(gv0, gvslice[:, 0]))
            # Check trait 1
            gv1 = gv1.reshape((len(w), 2))
            gv1 = np.sum(gv1, axis=1)
            self.assertTrue(np.allclose(gv1, gvslice[:, 1]))

            # Now, we check fitness
            for i, j, w in zip(gv0, gv1, md["w"]):
                d0 = np.power(i - self.zopt, 2.0)
                d1 = np.power(j - 0.0, 2.0)
                self.assertTrue(np.isclose(np.exp(-(d0 + d1) / (2.0 * self.VS)), w))

    def test_pickling(self):
        pp = pickle.dumps(self.pop, -1)
        up = pickle.loads(pp)
        gv = self.pop.genetic_values
        upgv = up.genetic_values
        self.assertTrue(np.array_equal(gv, upgv))
        gv = self.pop.ancient_sample_genetic_values
        upgv = up.ancient_sample_genetic_values
        self.assertTrue(np.array_equal(gv, upgv))

    def test_pickle_to_file(self):
        with tempfile.NamedTemporaryFile(delete=False) as tfile:
            self.pop.pickle_to_file(tfile)
            tfile.seek(0)
            pop2 = fwdpy11.DiploidPopulation.load_from_pickle_file(tfile)
        gv = self.pop.genetic_values
        gv2 = pop2.genetic_values
        self.assertTrue(np.array_equal(gv, gv2))
        gv = self.pop.ancient_sample_genetic_values
        gv2 = pop2.ancient_sample_genetic_values
        self.assertTrue(np.array_equal(gv, gv2))


class TestWithFirstGenerationPreserved(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.params, self.rng, self.pop = set_up_quant_trait_model()
        preserver = PreserveN(10 * self.pop.N, 10)
        fwdpy11.evolvets(
            self.rng,
            self.pop,
            self.params,
            100,
            preserver,
            record_gvalue_matrix=True,
            preserve_first_generation=True,
        )

    def test_alive_genetic_values(self):
        for i, j in zip(self.pop.genetic_values.flatten(), self.pop.diploid_metadata):
            self.assertEqual(i, j.g)

    def test_alive_genetic_value_reconstruction(self):
        ti = fwdpy11.TreeIterator(
            self.pop.tables, self.pop.alive_nodes, update_samples=True
        )
        gv = np.zeros(2 * self.pop.N)
        for t in ti:
            for m in t.mutations():
                for b in t.samples_below(m.node):
                    gv[b] += self.pop.mutations[m.key].s
        gv = gv.reshape((self.pop.N, 2))
        gv = np.sum(gv, axis=1)
        md = np.array(self.pop.diploid_metadata, copy=False)
        self.assertTrue(np.allclose(gv, md["g"]))

    def test_ancient_sample_genetic_values(self):
        for i, j in zip(
            self.pop.ancient_sample_genetic_values, self.pop.ancient_sample_metadata
        ):
            self.assertEqual(i, j.g)

    def test_ancient_sample_genetic_value_reconstruction(self):
        as_gv = self.pop.ancient_sample_genetic_values
        amd = np.array(self.pop.ancient_sample_metadata, copy=False)
        nt = np.array(self.pop.tables.nodes, copy=False)
        ancient_nodes = amd["nodes"][:, 0]
        ancient_node_times = nt["time"][ancient_nodes]
        for time, n, md in self.pop.sample_timepoints(include_alive=False):
            w = np.where(ancient_node_times == time)[0]
            gvslice = as_gv[w].flatten()
            ti = fwdpy11.TreeIterator(self.pop.tables, n, update_samples=True)
            gv = np.zeros(len(n))
            node_map = np.array([np.iinfo(np.int32).max] * len(nt), dtype=np.int32)
            for i, j in enumerate(n):
                node_map[j] = i
            for t in ti:
                for m in t.mutations():
                    for b in t.samples_below(m.node):
                        gv[node_map[b]] += self.pop.mutations[m.key].s
            gv = gv.reshape((len(w), 2))
            gv = np.sum(gv, axis=1)
            self.assertTrue(np.allclose(gv, gvslice))

    def test_pickling(self):
        pp = pickle.dumps(self.pop, -1)
        up = pickle.loads(pp)
        gv = self.pop.genetic_values
        upgv = up.genetic_values
        self.assertTrue(np.array_equal(gv, upgv))
        gv = self.pop.ancient_sample_genetic_values
        upgv = up.ancient_sample_genetic_values
        self.assertTrue(np.array_equal(gv, upgv))

    def test_pickle_to_file(self):
        with tempfile.NamedTemporaryFile(delete=False) as tfile:
            self.pop.pickle_to_file(tfile)
            tfile.seek(0)
            pop2 = fwdpy11.DiploidPopulation.load_from_pickle_file(tfile)
        gv = self.pop.genetic_values
        gv2 = pop2.genetic_values
        self.assertTrue(np.array_equal(gv, gv2))
        gv = self.pop.ancient_sample_genetic_values
        gv2 = pop2.ancient_sample_genetic_values
        self.assertTrue(np.array_equal(gv, gv2))


if __name__ == "__main__":
    unittest.main()
