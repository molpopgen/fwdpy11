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
import pickle
import unittest

import numpy as np

import fwdpy11


class TestPickling(unittest.TestCase):
    """
    Test pickling a couple of models that use all of the features.

    Note: test is limited as not all parameter types support __eq__ (yet)
    """

    @classmethod
    def setUp(self):
        self.pdict = {
            "nregions": [],
            "sregions": [fwdpy11.ExpS(0, 1, 1, -0.1, 0.25)],
            "recregions": [fwdpy11.PoissonInterval(0, 1, 1e-2)],
            "rates": [0, 0.01, None],
            "demography": fwdpy11.DiscreteDemography(),
            "simlen": 10,
            "gvalue": fwdpy11.Multiplicative(2.0),
        }

    def test_pickle_IM_with_selfing(self):
        import fwdpy11.demographic_models.IM as IM

        temp = IM.two_deme_IM(1000, 100, 0.25, [0.5, 3], [0.1, 0.25], as_dict=True)
        temp[0]["set_selfing_rates"] = [fwdpy11.SetSelfingRate(0, 0, 0.1)]
        im = fwdpy11.DiscreteDemography(**temp[0])
        self.pdict["demography"] = im
        params = fwdpy11.ModelParams(**self.pdict)
        pi = pickle.dumps(params)
        up = pickle.loads(pi)
        self.assertEqual(params.demography, up.demography)

    def test_pickle_IM_with_whole_migmatrix_reset(self):
        import fwdpy11.demographic_models.IM as IM

        temp = IM.two_deme_IM(1000, 100, 0.25, [0.5, 3], [0.1, 0.25], as_dict=True)
        temp[0]["set_migration_rates"].append(
            fwdpy11.SetMigrationRates(10000 + 10, None, np.eye(2))
        )
        im = fwdpy11.DiscreteDemography(**temp[0])
        self.pdict["demography"] = im
        params = fwdpy11.ModelParams(**self.pdict)
        pi = pickle.dumps(params)
        up = pickle.loads(pi)
        self.assertEqual(params.demography, up.demography)

    def test_pickle_tennessen(self):
        import fwdpy11.demographic_models.human as human

        tennessen = human.tennessen()
        self.pdict["demography"] = tennessen[0]
        params = fwdpy11.ModelParams(**self.pdict)
        pi = pickle.dumps(params)
        up = pickle.loads(pi)
        self.assertEqual(params.demography, up.demography)

    def test_pickle_model_with_scaled_migration_matrix(self):
        m = fwdpy11.DiscreteDemography(migmatrix=(np.eye(2), True))
        self.pdict["demography"] = m
        params = fwdpy11.ModelParams(**self.pdict)
        pi = pickle.dumps(params)
        up = pickle.loads(pi)
        self.assertEqual(params.demography, up.demography)


if __name__ == "__main__":
    unittest.main()
