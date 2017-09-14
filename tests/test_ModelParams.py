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

import fwdpy11 as fp11
import fwdpy11.model_params as fp11mp
import fwdpy11.multilocus as ml
import fwdpy11.wright_fisher_qtrait as wfq
from fwdpy11.fitness import SlocusAdditive
from quick_pops import quick_slocus_qtrait_pop_params
import numpy as np
import pickle


class testSlocusParams(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        # We do not assign certain attributes,
        # which allows us to test the defaults.
        self.pdict = {'nregions': [],
                      'sregions': [],
                      'recregions': [],
                      'rates': [0, 0, 0],
                      'demography': np.array([100] * 100, dtype=np.uint32)
                      }

    def test_SlocusParamsDefaults(self):
        from fwdpy11.fitness import SlocusMult
        m = fp11mp.SlocusParams(**self.pdict)
        self.assertEqual(type(m.gvalue), SlocusMult)

    def test_SlocusParamsQDefaults(self):
        self.pdict['prune_selected'] = False  # Suppresses a warning
        from fwdpy11.trait_values import SlocusAdditiveTrait
        from fwdpy11.wright_fisher_qtrait import GSS
        m = fp11mp.SlocusParamsQ(**self.pdict)
        self.assertEqual(type(m.gvalue), SlocusAdditiveTrait)
        self.assertEqual(type(m.trait2w), GSS)
        self.assertEqual(m.trait2w.VS, 1.0)
        self.assertEqual(m.trait2w.O, 0.0)

    def test_set_nregions(self):
        self.setUpClass()
        m = fp11mp.SlocusParams(**self.pdict)
        try:
            m.nregions = [fp11.Region(0, 1, 1)]
        except:
            self.fail("unexpected exception")
        # An empty list is ok
        try:
            m.nregions = []
        except:
            self.fail("unexpected exception")
        # TypeError b/c not iterable:
        with self.assertRaises(TypeError):
            m.nregions = fp11.Region(0, 1, 1)
        with self.assertRaises(TypeError):
            m.nregions = [fp11.ExpS(0, 1, 1, 0.25)]

    def test_set_recregions(self):
        self.setUpClass()
        m = fp11mp.SlocusParams(**self.pdict)
        try:
            m.recregions = [fp11.Region(0, 1, 1)]
        except:
            self.fail("unexpected exception")
        # An empty list is ok
        try:
            m.recregions = []
        except:
            self.fail("unexpected exception")
        # TypeError b/c not iterable:
        with self.assertRaises(TypeError):
            m.recregions = fp11.Region(0, 1, 1)
        with self.assertRaises(TypeError):
            m.recregions = [fp11.ExpS(0, 1, 1, 0.25)]

    def test_set_sregions(self):
        self.setUpClass()
        m = fp11mp.SlocusParams(**self.pdict)
        try:
            m.sregions = [fp11.GaussianS(0, 1, 1, 0.1)]
        except:
            self.fail("unexpected exception")
        # An empty list is ok
        try:
            m.sregions = []
        except:
            self.fail("unexpected exception")
        with self.assertRaises(TypeError):
            m.sregions = fp11.GaussianS(0, 1, 1, 0.2)
        with self.assertRaises(TypeError):
            m.sregions = [fp11.Region(0, 1, 1)]

    def test_init_demog(self):
        self.setUpClass()
        try:
            m = fp11mp.SlocusParams(**self.pdict)
            m
        except:
            self.fail("unexpected exception")

    def test_init_demog_list(self):
        self.pdict['demography'] = [100] * 10
        with self.assertRaises(ValueError):
            m = fp11mp.SlocusParams(**self.pdict)
            m

    def test_init_bad_popsizes(self):
        self.pdict['demography'] = np.array([100] * 10 + [-1])
        with self.assertRaises(ValueError):
            m = fp11mp.SlocusParams(**self.pdict)
            m

    def test_bad_rates(self):
        self.setUpClass()
        self.pdict['mutrate_n'] = None
        with self.assertRaises(ValueError):
            m = fp11mp.SlocusParams(**self.pdict)
            m
        self.pdict['mutrate_n'] = -0.01
        with self.assertRaises(ValueError):
            m = fp11mp.SlocusParams(**self.pdict)
            m

    def test_warning_for_zero_rate(self):
        # We pass in some sregions but leave mutrate_s == 0.
        self.setUpClass()
        self.pdict['sregions'] = [fp11.GammaS(0, 1, 1, -0.1, 0.10, 1.0)]
        with self.assertWarns(UserWarning):
            m = fp11mp.SlocusParams(**self.pdict)
            m.validate()

    def test_set_rates_list(self):
        self.setUpClass()
        self.pdict['rates'] = [0., 1e-3, 0.25]
        try:
            m = fp11mp.SlocusParams(**self.pdict)
            m
        except:
            self.fail("unexpected exception")

    def test_set_rates_tuple(self):
        m = fp11mp.SlocusParams()
        m.rates = (0., 1e-3, 0.25)
        self.assertEqual(m.mutrate_n, 0.0)
        self.assertEqual(m.mutrate_s, 1e-3)
        self.assertEqual(m.recrate, 0.25)
        with self.assertRaises(ValueError):
            m.rates = (0., -1e-3, 0.25)

    def test_set_rates_dict(self):
        m = fp11mp.SlocusParams()
        with self.assertRaises(ValueError):
            m.rates = {'mutrate_n': 1e-3}
        with self.assertRaises(ValueError):
            m.rates = {'mutrate_n': -1e-3,
                       'mutrate_s': 1e-3,
                       'recrate': 1e-3}
        with self.assertRaises(ValueError):
            m.rates = {'monkeys_n': 1e-3}

    def test_set_gvalue(self):
        m = fp11mp.SlocusParams()
        from fwdpy11.fitness import SlocusAdditive
        m.gvalue = SlocusAdditive()

    def test_set_bad_gvalue(self):
        m = fp11mp.SlocusParams()
        with self.assertRaises(ValueError):
            m.gvalue = []

    def test_invalid_kwargs(self):
        from fwdpy11.fitness import SlocusAdditive
        with self.assertRaises(ValueError):
            m = fp11mp.SlocusParams(fitness=SlocusAdditive())
            m


class testMlocusParams(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.rng = fp11.GSLrng(42)
        nregions = [[fp11.Region(0, 1, 1)], [fp11.Region(1, 2, 1)]]
        sregions = [
            [fp11.ExpS(0, 1, 1, -0.1)],
            [fp11.GammaS(1, 2, 1, -0.1, 0.5)]]
        recregions = nregions
        interlocus = ml.binomial_rec([0.5])
        region_rates = [1e-3, 1e-3]
        genetic_value = ml.MultiLocusGeneticValue([SlocusAdditive()] * 2)
        nlist = np.array([100] * 10, dtype=np.uint32)
        self.param_dict = {'nregions': nregions,
                           'sregions': sregions,
                           'recregions': recregions,
                           'interlocus': interlocus,
                           'mutrates_n': region_rates,
                           'mutrates_s': region_rates,
                           'recrates': region_rates,
                           'gvalue': genetic_value,
                           'demography': nlist}
        self.m = fp11mp.MlocusParams(**self.param_dict)

    def test_assign_gvalue(self):
        self.m.gvalue = ml.MultiLocusGeneticValue([SlocusAdditive()] * 2)

    def test_assign_gvalue_from_list(self):
        self.m.gvalue = [SlocusAdditive()] * 2

    def test_assign_gvalue_from_bad_input(self):
        with self.assertRaises(ValueError):
            self.m.gvalue = [2.0] * 2

    def test_assign_aggregator(self):
        self.m.aggregator = ml.AggAddTrait()

    def test_assign_aggregator_non_callable(self):
        with self.assertRaises(ValueError):
            self.m.aggregator = 2.0
    # Now, tests that check the validate(), too

    def test_construction(self):
        self.setUpClass()
        try:
            self.m.validate()
        except:
            self.fail("unexpected exception")

    def test_construction_unequal_region_lens_1(self):
        self.param_dict['nregions'] = [fp11.Region(0, 1, 1)]
        with self.assertRaises(TypeError):
            self.m = fp11mp.MlocusParams(**self.param_dict)

    def test_construction_unequal_region_lens_2(self):
        self.param_dict['sregions'] = [fp11.ExpS(0, 1, 1, 0.25)]
        with self.assertRaises(TypeError):
            self.m = fp11mp.MlocusParams(**self.param_dict)

    def test_construction_unequal_region_lens_3(self):
        self.param_dict['recregions'] = [fp11.Region(0, 1, 1)]
        with self.assertRaises(TypeError):
            self.m = fp11mp.MlocusParams(**self.param_dict)

    def test_empty_sregion_exception(self):
        self.setUpClass()
        self.param_dict['sregions'] = [[fp11.ExpS(0, 1, 1, 0.25)], []]
        self.m = fp11mp.MlocusParams(**self.param_dict)
        with self.assertRaises(ValueError):
            self.m.validate()

    def test_empty_nregion_exception(self):
        self.setUpClass()
        self.param_dict['nregions'] = [[], []]
        self.m = fp11mp.MlocusParams(**self.param_dict)
        with self.assertRaises(ValueError):
            self.m.validate()

    def test_empty_recregion_exception(self):
        self.setUpClass()
        self.param_dict['recregions'] = [[fp11.Region(0, 1, 1)], []]
        self.m = fp11mp.MlocusParams(**self.param_dict)
        with self.assertRaises(ValueError):
            self.m.validate()

    def test_pickle_pdict(self):
        self.setUpClass()
        p = pickle.dumps(self.param_dict)
        up = pickle.loads(p)
        up

    def test_pickle_params(self):
        self.setUpClass()
        self.m = fp11mp.MlocusParams(**self.param_dict)
        p = pickle.dumps(self.m)
        up = pickle.loads(p)
        up

    def test_deprecated_keyword(self):
        self.setUpClass()
        self.m = fp11mp.MlocusParams(**self.param_dict)
        with self.assertWarns(DeprecationWarning):
            self.m.agg = ml.AggAddTrait()


class testMlocusParamsDefaults(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        nregions = [[fp11.Region(0, 1, 1)], [fp11.Region(1, 2, 1)]]
        sregions = [
            [fp11.ExpS(0, 1, 1, -0.1)],
            [fp11.GammaS(1, 2, 1, -0.1, 0.5)]]
        recregions = nregions
        interlocus = ml.binomial_rec([0.5])
        region_rates = [1e-3, 1e-3]
        nlist = np.array([100] * 10, dtype=np.uint32)
        self.pdict = {'nregions': nregions,
                      'sregions': sregions,
                      'recregions': recregions,
                      'interlocus': interlocus,
                      'mutrates_n': region_rates,
                      'mutrates_s': region_rates,
                      'recrates': region_rates,
                      'demography': nlist
                      }

    def test_MlocusParamsDefaults(self):
        m = fp11mp.MlocusParams(**self.pdict)
        from fwdpy11.multilocus import AggMultFitness
        from fwdpy11.multilocus import MultiLocusGeneticValue
        from fwdpy11.fitness import SlocusMult
        self.assertEqual(type(m.aggregator), AggMultFitness)
        self.assertEqual(type(m.gvalue), MultiLocusGeneticValue)
        for i in m.gvalue.fitness_functions:
            self.assertEqual(type(i), SlocusMult)

    def test_MlocusParamsQDefaults(self):
        self.pdict['prune_selected'] = True
        m = fp11mp.MlocusParamsQ(**self.pdict)
        from fwdpy11.multilocus import AggMultTrait
        from fwdpy11.multilocus import MultiLocusGeneticValue
        from fwdpy11.trait_values import SlocusAdditiveTrait
        from fwdpy11.wright_fisher_qtrait import GSS
        self.assertEqual(type(m.aggregator), AggMultTrait)
        self.assertEqual(type(m.gvalue), MultiLocusGeneticValue)
        self.assertEqual(type(m.trait2w), GSS)
        self.assertEqual(m.trait2w.VS, 1.0)
        self.assertEqual(m.trait2w.O, 0.0)
        for i in m.gvalue.fitness_functions:
            self.assertEqual(type(i), SlocusAdditiveTrait)

    def test_MlocusParams_cannot_set_gvalue(self):
        self.pdict.pop('sregions')
        m = fp11mp.MlocusParams(**self.pdict)
        self.assertIs(m.gvalue, None)
        m = fp11mp.MlocusParamsQ(**self.pdict)
        self.assertIs(m.gvalue, None)


class test_SlocusParamsQEvolve(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.pop, self.params = quick_slocus_qtrait_pop_params()

    def test_evolve(self):
        rng = fp11.GSLrng(42)
        wfq.evolve(rng, self.pop, self.params)


if __name__ == "__main__":
    unittest.main()
