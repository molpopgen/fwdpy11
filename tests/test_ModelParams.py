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
from fwdpy11.fitness import SlocusAdditive
import numpy as np


class testModelParamsConstructor(unittest.TestCase):
    def test_init(self):
        m = fp11mp.ModelParams()
        try:
            m.validate()
        except ValueError:
            self.fail("ValueError encountered")

    def test_nregions(self):
        m = fp11mp.ModelParams(nregions=[fp11.Region(0, 1, 1)])
        self.assertIsInstance(m.nregions, list)
        self.assertIsInstance(m.nregions[0], fp11.Region)
        try:
            m.validate()
        except ValueError:
            self.fail("ValueError encountered")

    def test_recregions(self):
        m = fp11mp.ModelParams(recregions=[fp11.Region(0, 1, 1)])
        self.assertIsInstance(m.recregions, list)
        self.assertIsInstance(m.recregions[0], fp11.Region)
        try:
            m.validate()
        except ValueError:
            self.fail("ValueError encountered")

    def test_sregions(self):
        m = fp11mp.ModelParams(sregions=[fp11.ExpS(0, 1, 1, 0.25)])
        self.assertIsInstance(m.sregions, list)
        self.assertIsInstance(m.sregions[0], fp11.ExpS)
        try:
            m.validate()
        except ValueError:
            self.fail("ValueError encountered")


class testModelParamsSetting(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.m = fp11mp.ModelParams()

    def test_nregions(self):
        self.m.nregions = [fp11.Region(0, 1, 1)]
        self.assertIsInstance(self.m.nregions, list)
        self.assertIsInstance(self.m.nregions[0], fp11.Region)
        try:
            self.m.validate()
        except ValueError:
            self.fail("ValueError encountered")

    def test_recregions(self):
        self.m.recregions = [fp11.Region(0, 1, 1)]
        self.assertIsInstance(self.m.recregions, list)
        self.assertIsInstance(self.m.recregions[0], fp11.Region)
        try:
            self.m.validate()
        except ValueError:
            self.fail("ValueError encountered")

    def test_sregions(self):
        self.m.sregions = [fp11.ExpS(0, 1, 1, 0.25)]
        self.assertIsInstance(self.m.sregions, list)
        self.assertIsInstance(self.m.sregions[0], fp11.ExpS)
        try:
            self.m.validate()
        except ValueError:
            self.fail("ValueError encountered")


class testModelParamsBadInputData(unittest.TestCase):
    def test_not_list(self):
        m = fp11mp.ModelParams(nregions=fp11.Region(0, 1, 1))
        with self.assertRaises(ValueError):
            m.validate()


class testSlocusParams(unittest.TestCase):
    def test_init_demog(self):
        m = fp11mp.SlocusParams(
            demography=np.array([100] * 10, dtype=np.uint32))
        m.validate()

    def test_init_demog_list(self):
        with self.assertRaises(ValueError):
            m = fp11mp.SlocusParams(demography=[100] * 10)

    def test_init_bad_popsizes(self):
        with self.assertRaises(ValueError):
            m = fp11mp.SlocusParams(demography=np.array([100] * 10 + [-1]))

    def test_bad_rates(self):
        with self.assertRaises(ValueError):
            m = fp11mp.SlocusParams(mutrate_n=None)

    def test_warning_for_zero_rate(self):
        # We pass in some sregions but leave mutrate_s == 0
        m = fp11mp.SlocusParams(
            sregions=[fp11.GammaS(0, 1, 1, -0.1, 0.10, 1.0)])
        m.demography = np.array([100] * 100, dtype=np.uint32)
        with self.assertWarns(UserWarning):
            m.validate()

    def test_set_rates_list(self):
        m = fp11mp.SlocusParams()
        m.rates = [0., 1e-3, 0.25]
        with self.assertRaises(ValueError):
            m.rates = [0., -1e-3, 0.25]

    def test_set_rates_tuple(self):
        m = fp11mp.SlocusParams()
        m.rates = (0., 1e-3, 0.25)
        with self.assertRaises(ValueError):
            m.rates = (0., -1e-3, 0.25)

    def test_set_rates_dict(self):
        m = fp11mp.SlocusParams()
        m.rates = {'mutrate_n': 1e-3}
        with self.assertRaises(ValueError):
            m.rates = {'mutrate_n': -1e-3}
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


class testMlocusParams(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.m = fp11mp.MlocusParams()
        self.rng = fp11.GSLrng(42)
        nregions = [[fp11.Region(0, 1, 1)], [fp11.Region(1, 2, 1)]]
        sregions = [
            [fp11.ExpS(0, 1, 1, -0.1)],
            [fp11.GammaS(1, 2, 1, -0.1, 0.5)]]
        recregions = nregions
        interlocus = ml.binomial_rec(self.rng, [0.5])
        region_rates = [1e-3, 1e-3]
        genetic_value = ml.MultiLocusGeneticValue([SlocusAdditive()] * 2)
        self.param_dict = {'nregions': nregions,
                           'sregions': sregions,
                           'recregions': recregions,
                           'interlocus': interlocus,
                           'mutrates_n': region_rates,
                           'mutrates_s': region_rates,
                           'recrates': region_rates,
                           'gvalue': genetic_value}

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
        self.m = fp11mp.MlocusParams(**self.param_dict)
        try:
            self.m.validate()
        except:
            self.fail("unexpected exception")

    def test_construction_unequal_region_lens_1(self):
        self.param_dict['nregions'] = [fp11.Region(0, 1, 1)]
        self.m = fp11mp.MlocusParams(**self.param_dict)
        with self.assertRaises(ValueError):
            self.m.validate()

    def test_construction_unequal_region_lens_2(self):
        self.param_dict['sregions'] = [fp11.ExpS(0, 1, 1, 0.25)]
        self.m = fp11mp.MlocusParams(**self.param_dict)
        with self.assertRaises(ValueError):
            self.m.validate()

    def test_construction_unequal_region_lens_3(self):
        self.param_dict['recregions'] = [fp11.Region(0, 1, 1)]
        self.m = fp11mp.MlocusParams(**self.param_dict)
        with self.assertRaises(ValueError):
            self.m.validate()

if __name__ == "__main__":
    unittest.main()

