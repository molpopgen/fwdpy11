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

import fwdpy11
import unittest
import numpy as np


def setup_pdict(demog, simlen):
    nregions = []
    sregions = []
    recregions = []

    pdict = {'nregions': nregions,
             'sregions': sregions,
             'recregions': recregions,
             'rates': (0, 0, None),
             'gvalue': fwdpy11.Multiplicative(2.),
             'demography': demog,
             'simlen': simlen,
             'prune_selected': True
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
        self.d, self.t1, self.t2 = two_deme_IM(self.Nanc, 0.1, 0.7,
                                               (self.N0, self.N1),
                                               (1e-2, 0.25),
                                               burnin=1.0)
        self.pop, self.rng = setup_pop_rng(self.Nanc)

    def test_complete_sim(self):
        pdict = setup_pdict(self.d, self.t1+self.t2)
        params = fwdpy11.ModelParams(**pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 10)
        self.assertEqual(self.pop.generation, params.simlen)
        deme_sizes = self.pop.deme_sizes(as_dict=True)
        self.assertEqual(deme_sizes[0], np.rint(self.N0*self.Nanc).astype(int))
        self.assertEqual(deme_sizes[1], np.rint(self.N1*self.Nanc).astype(int))

    def test_evolve_in_two_steps(self):
        pdict = setup_pdict(self.d, self.t1)
        params = fwdpy11.ModelParams(**pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 10)
        self.assertEqual(self.pop.generation, self.t1)
        deme_sizes = self.pop.deme_sizes(as_dict=True)
        self.assertEqual(deme_sizes[0], self.Nanc)
        self.assertTrue(1 not in deme_sizes)

        params.simlen = self.t2
        fwdpy11.evolvets(self.rng, self.pop, params, 10)
        self.assertEqual(self.pop.generation, self.t1+self.t2)
        deme_sizes = self.pop.deme_sizes(as_dict=True)
        self.assertEqual(deme_sizes[0], np.rint(self.N0*self.Nanc).astype(int))
        self.assertEqual(deme_sizes[1], np.rint(self.N1*self.Nanc).astype(int))

    def test_evolve_in_two_steps_restart_with_two_demes(self):
        deltat = 2
        pdict = setup_pdict(self.d, self.t1+deltat)
        params = fwdpy11.ModelParams(**pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 10)
        self.assertEqual(self.pop.generation, self.t1+deltat)
        deme_sizes = self.pop.deme_sizes(as_dict=True)
        self.assertTrue(0 in deme_sizes)
        self.assertTrue(1 in deme_sizes)

        params.simlen = self.t2-deltat
        fwdpy11.evolvets(self.rng, self.pop, params, 10)
        self.assertEqual(self.pop.generation, self.t1+self.t2)
        deme_sizes = self.pop.deme_sizes(as_dict=True)
        self.assertEqual(deme_sizes[0], np.rint(self.N0*self.Nanc).astype(int))
        self.assertEqual(deme_sizes[1], np.rint(self.N1*self.Nanc).astype(int))


if __name__ == "__main__":
    unittest.main()
