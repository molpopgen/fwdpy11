#
# Copyright (C) 2020 Kevin Thornton <krthornt@uci.edu>
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

import scipy.stats

import fwdpy11
import sregion_cdf


class TestSregionFromCDF(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        # Test everything w.r.to a value of 0.1
        # from the underlying cdf
        self.P = 0.1

    def get_deviate(self, dist):
        return sregion_cdf.from_mvnorm(0.0, self.P, dist)

    def test_ExpS(self):
        mean = 0.2
        e = fwdpy11.ExpS(0, 1, 1, mean)
        d = self.get_deviate(e)
        exp = scipy.stats.expon(0, mean)
        self.assertAlmostEqual(self.P, exp.cdf(d))

    def test_ExpS_with_scaling(self):
        mean = 0.2
        scaling = 1e3
        e = fwdpy11.ExpS(0, 1, 1, mean * scaling, scaling=scaling)
        d = self.get_deviate(e)
        exp = scipy.stats.expon(0, mean)
        self.assertAlmostEqual(self.P, exp.cdf(d))

    def test_GaussianS(self):
        sd = 0.333
        g = fwdpy11.GaussianS(0, 1, 1, sd)
        d = self.get_deviate(g)
        norm = scipy.stats.norm(0, sd)
        self.assertAlmostEqual(self.P, norm.cdf(d))

    def test_GaussianS_with_scaling(self):
        sd = 0.333
        scaling = 1e3
        g = fwdpy11.GaussianS(0, 1, 1, sd * scaling, scaling=scaling)
        d = self.get_deviate(g)
        norm = scipy.stats.norm(0, sd)
        self.assertAlmostEqual(self.P, norm.cdf(d))

    def test_GammaS(self):
        # a = shape
        # b = mean/shape
        mean = 5
        shape = 2.6
        g = fwdpy11.GammaS(0, 1, 1, mean=mean, shape_parameter=shape)
        d = self.get_deviate(g)
        gamma = scipy.stats.gamma(a=shape, scale=mean / shape)
        self.assertAlmostEqual(self.P, gamma.cdf(d))

    def test_GammaS_with_scaling(self):
        # a = shape
        # b = mean/shape
        scaling = 1e3
        mean = 5 * scaling
        shape = 2.6
        g = fwdpy11.GammaS(0, 1, 1, mean=mean, shape_parameter=shape, scaling=scaling)
        d = self.get_deviate(g)
        gamma = scipy.stats.gamma(a=shape, scale=mean / (scaling * shape))
        self.assertAlmostEqual(self.P, gamma.cdf(d))

    def test_UniformS(self):
        lo = 0.1
        hi = 0.47
        u = fwdpy11.UniformS(0, 1, 1, lo, hi)
        d = self.get_deviate(u)
        uniform = scipy.stats.uniform(loc=lo, scale=hi - lo)
        self.assertAlmostEqual(self.P, uniform.cdf(d))

    def test_UniformS_with_scaling(self):
        lo = 100
        hi = 200
        scaling = 1e3
        u = fwdpy11.UniformS(0, 1, 1, lo, hi, scaling=scaling)
        d = self.get_deviate(u)
        uniform = scipy.stats.uniform(loc=lo, scale=hi - lo)
        self.assertAlmostEqual(self.P, uniform.cdf(d * scaling))


if __name__ == "__main__":
    unittest.main()
