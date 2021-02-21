#
# Copyright (C) 2021 Kevin Thornton <krthornt@uci.edu>
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

import fwdpy11
import numpy as np
import pytest


@pytest.mark.parametrize(
    "t", [(fwdpy11.FixedDominance, "h"), (fwdpy11.ExponentialDominance, "m")]
)
@pytest.mark.parametrize("x", [1.0, -0.25])
def test_single_param_construction(t, x):
    _ = t[0](x)
    kwargs = {t[1]: x}
    _ = t[0](**kwargs)


@pytest.mark.parametrize(
    "t", [(fwdpy11.FixedDominance, "h"), (fwdpy11.ExponentialDominance, "m")]
)
@pytest.mark.parametrize("x", [np.nan, np.inf])
def test_single_param_bad_construction(t, x):
    with pytest.raises(ValueError):
        _ = t[0](x)
    with pytest.raises(ValueError):
        kwargs = {t[1]: x}
        _ = t[0](**kwargs)


@pytest.mark.parametrize("x", [(0.0, 1.0), (-0.25, 1.0), (-1.0, 0.0), (-1.0, -1 / 3)])
def test_UniformDominance_construction(x):
    _ = fwdpy11.UniformDominance(lo=x[0], hi=x[1])
    _ = fwdpy11.UniformDominance(x[0], x[1])


@pytest.mark.parametrize(
    "x", [(np.inf, 1.0), (np.nan, 1.0), (1.0, np.inf), (np.nan, 1.0), (1.0, 0.0)]
)
def test_UniformDominance_bad_construction(x):
    with pytest.raises(ValueError):
        _ = fwdpy11.UniformDominance(lo=x[0], hi=x[1])


@pytest.mark.parametrize("k", [1e-3, 1.0, 10.0])
def test_LargeEffectExponentiallyRecessive_construction(k):
    _ = fwdpy11.LargeEffectExponentiallyRecessive(k)
    _ = fwdpy11.LargeEffectExponentiallyRecessive(k=k)


@pytest.mark.parametrize("k", [-1e-3, np.nan, np.inf])
def test_LargeEffectExponentiallyRecessive_bad_construction(k):
    with pytest.raises(ValueError):
        _ = fwdpy11.LargeEffectExponentiallyRecessive(k)
    with pytest.raises(ValueError):
        _ = fwdpy11.LargeEffectExponentiallyRecessive(k=k)


@pytest.mark.parametrize(
    "x",
    [
        fwdpy11.FixedDominance(h=0.667),
        fwdpy11.ExponentialDominance(m=-0.23),
        fwdpy11.LargeEffectExponentiallyRecessive(k=1.0),
        fwdpy11.LargeEffectExponentiallyRecessive(k=1.0, scaling=1.),
        fwdpy11.UniformDominance(lo=-2 / 3, hi=2 / 3),
    ],
)
def test_pickling(x):
    p = pickle.dumps(x, -1)
    up = pickle.loads(p)
    assert up == x


@pytest.mark.parametrize(
    "des",
    [
        (fwdpy11.ExpS, {"mean": -0.2}),
        (fwdpy11.ConstantS, {"s": -0.2}),
        (fwdpy11.GammaS, {"mean": 1, "shape_parameter": 1}),
        (fwdpy11.GaussianS, {"sd": 1}),
        (fwdpy11.LogNormalS, {"zeta": 1, "sigma": 1}),
        (fwdpy11.UniformS, {"lo": -0.1, "hi": -1e-3}),
        (fwdpy11.MultivariateGaussianEffects, {"cov_matrix": np.identity(2)}),
    ],
)
@pytest.mark.parametrize(
    "h",
    [
        fwdpy11.FixedDominance(h=0.667),
        fwdpy11.ExponentialDominance(m=-0.23),
        fwdpy11.LargeEffectExponentiallyRecessive(k=1.0),
        fwdpy11.LargeEffectExponentiallyRecessive(k=1.0, scaling=0.5),
        fwdpy11.UniformDominance(lo=-2 / 3, hi=2 / 3),
    ],
)
def test_build_DES_objects(des, h):
    x = des[0](beg=0, end=1, weight=1, h=h, **des[1])
    p = pickle.dumps(x)
    up = pickle.loads(p)
    assert x == up
