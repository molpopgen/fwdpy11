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

# TODO reimplement to test new API

# import fwdpy11 as fp11
# import fwdpy11 as fp11mp
# import fwdpy11.multilocus as ml
# import fwdpy11.wright_fisher_qtrait as wfq
# from fwdpy11.fitness import DiploidAdditive
# from quick_pops import quick_slocus_qtrait_pop_params
# import numpy as np
# import pickle
# 
# 
# class testDiploidParams(unittest.TestCase):
#     @classmethod
#     def setUpClass(self):
#         # We do not assign certain attributes,
#         # which allows us to test the defaults.
#         self.pdict = {'nregions': [],
#                       'sregions': [],
#                       'recregions': [],
#                       'rates': [0, 0, 0],
#                       'demography': np.array([100] * 100, dtype=np.uint32)
#                       }
# 
#     def test_DiploidParamsDefaults(self):
#         from fwdpy11.fitness import DiploidMult
#         m = fp11mp.DiploidParams(**self.pdict)
#         self.assertEqual(type(m.gvalue), DiploidMult)
# 
#     def test_DiploidParamsQDefaults(self):
#         self.pdict['prune_selected'] = False  # Suppresses a warning
#         from fwdpy11.trait_values import DiploidAdditiveTrait
#         from fwdpy11.wright_fisher_qtrait import GSS
#         m = fp11mp.DiploidParamsQ(**self.pdict)
#         self.assertEqual(type(m.gvalue), DiploidAdditiveTrait)
#         self.assertEqual(type(m.trait2w), GSS)
#         self.assertEqual(m.trait2w.VS, 1.0)
#         self.assertEqual(m.trait2w.O, 0.0)
# 
#     def test_set_nregions(self):
#         self.setUpClass()
#         m = fp11mp.DiploidParams(**self.pdict)
#         try:
#             m.nregions = [fp11.Region(0, 1, 1)]
#         except:
#             self.fail("unexpected exception")
#         # An empty list is ok
#         try:
#             m.nregions = []
#         except:
#             self.fail("unexpected exception")
#         # TypeError b/c not iterable:
#         with self.assertRaises(TypeError):
#             m.nregions = fp11.Region(0, 1, 1)
#         with self.assertRaises(TypeError):
#             m.nregions = [fp11.ExpS(0, 1, 1, 0.25)]
# 
#     def test_set_recregions(self):
#         self.setUpClass()
#         m = fp11mp.DiploidParams(**self.pdict)
#         try:
#             m.recregions = [fp11.Region(0, 1, 1)]
#         except:
#             self.fail("unexpected exception")
#         # An empty list is ok
#         try:
#             m.recregions = []
#         except:
#             self.fail("unexpected exception")
#         # TypeError b/c not iterable:
#         with self.assertRaises(TypeError):
#             m.recregions = fp11.Region(0, 1, 1)
#         with self.assertRaises(TypeError):
#             m.recregions = [fp11.ExpS(0, 1, 1, 0.25)]
# 
#     def test_set_sregions(self):
#         self.setUpClass()
#         m = fp11mp.DiploidParams(**self.pdict)
#         try:
#             m.sregions = [fp11.GaussianS(0, 1, 1, 0.1)]
#         except:
#             self.fail("unexpected exception")
#         # An empty list is ok
#         try:
#             m.sregions = []
#         except:
#             self.fail("unexpected exception")
#         with self.assertRaises(TypeError):
#             m.sregions = fp11.GaussianS(0, 1, 1, 0.2)
#         with self.assertRaises(TypeError):
#             m.sregions = [fp11.Region(0, 1, 1)]
# 
#     def test_init_demog(self):
#         self.setUpClass()
#         try:
#             m = fp11mp.DiploidParams(**self.pdict)
#             m
#         except:
#             self.fail("unexpected exception")
# 
#     def test_init_demog_list(self):
#         self.pdict['demography'] = [100] * 10
#         with self.assertRaises(ValueError):
#             m = fp11mp.DiploidParams(**self.pdict)
#             m
# 
#     def test_init_bad_popsizes(self):
#         self.pdict['demography'] = np.array([100] * 10 + [-1])
#         with self.assertRaises(ValueError):
#             m = fp11mp.DiploidParams(**self.pdict)
#             m
# 
#     def test_bad_rates(self):
#         self.setUpClass()
#         self.pdict['mutrate_n'] = None
#         with self.assertRaises(ValueError):
#             m = fp11mp.DiploidParams(**self.pdict)
#             m
#         self.pdict['mutrate_n'] = -0.01
#         with self.assertRaises(ValueError):
#             m = fp11mp.DiploidParams(**self.pdict)
#             m
# 
#     def test_warning_for_zero_rate(self):
#         # We pass in some sregions but leave mutrate_s == 0.
#         self.setUpClass()
#         self.pdict['sregions'] = [fp11.GammaS(0, 1, 1, -0.1, 0.10, 1.0)]
#         with self.assertWarns(UserWarning):
#             m = fp11mp.DiploidParams(**self.pdict)
#             m.validate()
# 
#     def test_set_rates_list(self):
#         self.setUpClass()
#         self.pdict['rates'] = [0., 1e-3, 0.25]
#         try:
#             m = fp11mp.DiploidParams(**self.pdict)
#             m
#         except:
#             self.fail("unexpected exception")
# 
#     def test_set_rates_tuple(self):
#         m = fp11mp.DiploidParams()
#         m.rates = (0., 1e-3, 0.25)
#         self.assertEqual(m.mutrate_n, 0.0)
#         self.assertEqual(m.mutrate_s, 1e-3)
#         self.assertEqual(m.recrate, 0.25)
#         with self.assertRaises(ValueError):
#             m.rates = (0., -1e-3, 0.25)
# 
#     def test_set_rates_dict(self):
#         m = fp11mp.DiploidParams()
#         with self.assertRaises(ValueError):
#             m.rates = {'mutrate_n': 1e-3}
#         with self.assertRaises(ValueError):
#             m.rates = {'mutrate_n': -1e-3,
#                        'mutrate_s': 1e-3,
#                        'recrate': 1e-3}
#         with self.assertRaises(ValueError):
#             m.rates = {'monkeys_n': 1e-3}
# 
#     def test_set_gvalue(self):
#         m = fp11mp.DiploidParams()
#         from fwdpy11.fitness import DiploidAdditive
#         m.gvalue = DiploidAdditive()
# 
#     def test_set_bad_gvalue(self):
#         m = fp11mp.DiploidParams()
#         with self.assertRaises(ValueError):
#             m.gvalue = []
# 
#     def test_invalid_kwargs(self):
#         from fwdpy11.fitness import DiploidAdditive
#         with self.assertRaises(ValueError):
#             m = fp11mp.DiploidParams(fitness=DiploidAdditive())
#             m
#
# if __name__ == "__main__":
#     pass
