#
# Copyright (C) 2019 Kevin Thornton <krthornt@uci.edu>
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

import os
import unittest

import numpy as np

import fwdpy11


def setup_pdict(demog, simlen):
    nregions = []
    sregions = []
    recregions = []

    pdict = {
        "nregions": nregions,
        "sregions": sregions,
        "recregions": recregions,
        "rates": (0, 0, None),
        "gvalue": fwdpy11.Multiplicative(2.0),
        "demography": demog,
        "simlen": simlen,
        "prune_selected": True,
    }
    return pdict


def setup_pop_rng(N=100, seed=42):
    pop = fwdpy11.DiploidPopulation(N, 1.0)
    rng = fwdpy11.GSLrng(seed)
    return pop, rng


class TestTwoDemeIMModel(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.Nanc = 100
        self.N0, self.N1 = 1.1, 3.0
        from fwdpy11.demographic_models.IM import two_deme_IM

        self.model = two_deme_IM(
            self.Nanc, 0.1, 0.7, (self.N0, self.N1), (1e-2, 0.25), burnin=1.0
        )
        self.pop, self.rng = setup_pop_rng(self.Nanc)

    def test_complete_sim(self):
        pdict = setup_pdict(
            self.model.model,
            self.model.metadata.split_time + self.model.metadata.gens_post_split,
        )
        params = fwdpy11.ModelParams(**pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 10)
        self.assertEqual(self.pop.generation, params.simlen)
        deme_sizes = self.pop.deme_sizes(as_dict=True)
        self.assertEqual(deme_sizes[0], np.rint(self.N0 * self.Nanc).astype(int))
        self.assertEqual(deme_sizes[1], np.rint(self.N1 * self.Nanc).astype(int))

    def test_evolve_in_two_steps(self):
        pdict = setup_pdict(self.model.model, self.model.metadata.split_time)
        params = fwdpy11.ModelParams(**pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 10)
        self.assertEqual(self.pop.generation, self.model.metadata.split_time)
        deme_sizes = self.pop.deme_sizes(as_dict=True)
        self.assertEqual(deme_sizes[0], int((1.0 - 0.7) * self.Nanc))

        pdict["simlen"] = self.model.metadata.gens_post_split
        params = fwdpy11.ModelParams(**pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 10)
        self.assertEqual(
            self.pop.generation,
            self.model.metadata.split_time + self.model.metadata.gens_post_split,
        )
        deme_sizes = self.pop.deme_sizes(as_dict=True)
        self.assertEqual(deme_sizes[0], np.rint(self.N0 * self.Nanc).astype(int))
        self.assertEqual(deme_sizes[1], np.rint(self.N1 * self.Nanc).astype(int))

    def test_evolve_in_two_steps_restart_with_two_demes(self):
        deltat = 2
        pdict = setup_pdict(self.model.model, self.model.metadata.split_time + deltat)
        params = fwdpy11.ModelParams(**pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 10)
        self.assertEqual(self.pop.generation, self.model.metadata.split_time + deltat)
        deme_sizes = self.pop.deme_sizes(as_dict=True)
        self.assertTrue(0 in deme_sizes)
        self.assertTrue(1 in deme_sizes)

        pdict["simlen"] = self.model.metadata.gens_post_split - deltat
        params = fwdpy11.ModelParams(**pdict)
        fwdpy11.evolvets(
            self.rng, self.pop, params, 10, check_demographic_event_timings=False
        )
        self.assertEqual(
            self.pop.generation,
            self.model.metadata.split_time + self.model.metadata.gens_post_split,
        )
        deme_sizes = self.pop.deme_sizes(as_dict=True)
        self.assertEqual(deme_sizes[0], np.rint(self.N0 * self.Nanc).astype(int))
        self.assertEqual(deme_sizes[1], np.rint(self.N1 * self.Nanc).astype(int))


class TestTwoDemeIMModelVariousInitMethods(unittest.TestCase):
    def test_init_migrates_tuple(self):
        from fwdpy11.demographic_models.IM import two_deme_IM

        model = two_deme_IM(1000, 0.1, 0.7, (1.1, 2.7), (1e-2, 0.25), burnin=1.0)
        self.assertEqual(model, model)

    def test_init_migrates_list(self):
        from fwdpy11.demographic_models.IM import two_deme_IM

        model = two_deme_IM(1000, 0.1, 0.7, (1.1, 2.7), [1e-2, 0.25], burnin=1.0)
        self.assertEqual(model, model)

    def test_init_migrates_numpy(self):
        from fwdpy11.demographic_models.IM import two_deme_IM

        model = two_deme_IM(
            1000, 0.1, 0.7, (1.1, 2.7), np.array([1e-2, 0.25]), burnin=1.0
        )
        self.assertEqual(model, model)


@unittest.skipIf(
    "RUN_EXPENSIVE_TESTS" not in os.environ or os.environ["RUN_EXPENSIVE_TESTS"] != "1",
    "Expensive test.",
)
class TestTennessenModel(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        from fwdpy11.demographic_models.human import tennessen

        self.demog = tennessen(0)
        self.simlen = self.demog.metadata["simlen"]
        self.Nref = self.demog.metadata["Nref"]
        self.pop = fwdpy11.DiploidPopulation(self.Nref, 1.0)
        self.pdict = setup_pdict(self.demog, self.simlen)
        self.params = fwdpy11.ModelParams(**self.pdict)
        self.rng = fwdpy11.GSLrng(915153)
        fwdpy11.evolvets(self.rng, self.pop, self.params, 100)

    def test_generation(self):
        self.assertEqual(self.pop.generation, self.simlen)

    def test_final_deme_sizes(self):
        md = np.array(self.pop.diploid_metadata, copy=False)
        deme_sizes = np.unique(md["deme"], return_counts=True)
        expected_deme_sizes = [420000, 512000]
        for i, j in zip(deme_sizes[1], expected_deme_sizes):
            self.assertEqual(i, j)


@unittest.skipIf(
    "RUN_EXPENSIVE_TESTS" not in os.environ or os.environ["RUN_EXPENSIVE_TESTS"] != "1",
    "Expensive test.",
)
class TestTennessenModelV1(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        from fwdpy11.demographic_models.human import tennessen

        self.demog = tennessen(0, fwdpy11.demographic_models.human.TennessenModel.V1)
        self.simlen = self.demog.metadata["simlen"]
        self.Nref = self.demog.metadata["Nref"]
        self.pop = fwdpy11.DiploidPopulation(self.Nref, 1.0)
        self.pdict = setup_pdict(self.demog, self.simlen)
        self.params = fwdpy11.ModelParams(**self.pdict)
        self.rng = fwdpy11.GSLrng(915153)
        fwdpy11.evolvets(self.rng, self.pop, self.params, 100)

    def test_generation(self):
        self.assertEqual(self.pop.generation, self.simlen)

    def test_final_deme_sizes(self):
        md = np.array(self.pop.diploid_metadata, copy=False)
        deme_sizes = np.unique(md["deme"], return_counts=True)
        expected_deme_sizes = [423125, 501425]
        for i, j in zip(deme_sizes[1], expected_deme_sizes):
            self.assertEqual(i, j)


@unittest.skipIf(
    "RUN_EXPENSIVE_TESTS" not in os.environ or os.environ["RUN_EXPENSIVE_TESTS"] != "1",
    "Expensive test.",
)
class TestJouganousThreeDemeModel(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        from fwdpy11.demographic_models.human import jouganous_three_deme

        self.demog = jouganous_three_deme(0)
        self.simlen = self.demog.metadata["simlen"]
        self.Nref = self.demog.metadata["Nref"]
        self.pop = fwdpy11.DiploidPopulation(self.Nref, 1.0)
        self.pdict = setup_pdict(self.demog, self.simlen)
        self.params = fwdpy11.ModelParams(**self.pdict)
        self.rng = fwdpy11.GSLrng(915153)
        fwdpy11.evolvets(self.rng, self.pop, self.params, 100)

    def test_generation(self):
        self.assertEqual(self.pop.generation, self.simlen)

    def test_final_deme_sizes(self):
        md = np.array(self.pop.diploid_metadata, copy=False)
        deme_sizes = np.unique(md["deme"], return_counts=True)
        expected_deme_sizes = [23721, 39611, 83681]
        for i, j in zip(deme_sizes[1], expected_deme_sizes):
            self.assertEqual(i, j)


if __name__ == "__main__":
    unittest.main()
