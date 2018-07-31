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

import fwdpy11
import fwdpy11.model_params
import fwdpy11.genetic_values
import fwdpy11.genetic_value_noise
import numpy as np
import pickle


class testGeneticValuesWithoutNoise(unittest.TestCase):
    @classmethod
    def setUp(self):
        """
        Set up a minimal fixture
        """
        self.pdict = {'nregions': [fwdpy11.Region(0, 1, 1)],
                      'sregions': [fwdpy11.ExpS(0, 1, 1, -1.0, 0.5)],
                      'recregions': [fwdpy11.Region(0, 1, 1)],
                      'rates': (1e-3, 1e-3, 1e-3),
                      'demography': np.array([100]*10, dtype=np.uint32),
                      }

    def testGeneticValueWithoutNoise(self):
        """
        Genetic value parameters must be stored as objects
        useable to compose an object instance.
        """
        self.pdict['gvalue'] = (fwdpy11.genetic_values.SlocusMult, (1,))
        p = fwdpy11.model_params.ModelParams(**self.pdict)
        x = p.make_gvalue()
        self.assertEqual(x.scaling, 1.0)

    def testGeneticValueWithoutNoiseException(self):
        """
        In response to GH issue #129
        """
        # This is a genetic value type of length 1.
        # Gvalue now has an odd type, etc.,
        # and accommodating this adds a lot more if/else
        # to ModelParams.make_gvalue, which will probably
        # already get longer in order to support dicts
        # (see next test).
        self.pdict['gvalue'] = fwdpy11.genetic_values.SlocusGBR
        p = fwdpy11.model_params.ModelParams(**self.pdict)
        with self.assertRaise(TypeError):
            """
            Type error b/c not of length 1
            """
            x = p.make_gvalue()

    def testGeneticValueWithoutNoiseDict(self):
        """
        Ideally, we would support
        construction via a dict, too, but
        the *tuple vs **dict is a bit annoying
        """
        self.pdict['gvalue'] = (
            fwdpy11.genetic_values.SlocusMult, {'scaling': 1.0})
        p = fwdpy11.model_params.ModelParams(**self.pdict)
        with self.assertRaises(TypeError):
            x = p.make_gvalue()
            self.assertEqual(x.scaling, 1.0)


if __name__ == "__main__":
    unittest.main()
