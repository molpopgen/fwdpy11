import pickle

import pytest

import fwdpy11


@pytest.mark.parametrize("discrete", [True, False])
def test_poisson_interval(discrete):
    pi = fwdpy11.PoissonInterval(0, 2, 1e-3, discrete=discrete)
    p = pickle.dumps(pi)
    up = pickle.loads(p)
    assert up.beg == pi.beg
    assert up.end == pi.end
    assert up.mean == pi.mean
    assert up.discrete == pi.discrete
    assert up.discrete == discrete

    if discrete is True:
        with pytest.raises(ValueError):
            fwdpy11.PoissonInterval(0, 1, 1e-3, discrete=discrete)
        with pytest.raises(TypeError):
            fwdpy11.PoissonInterval(0.0, 2, 1e-3, discrete=discrete)
        with pytest.raises(TypeError):
            fwdpy11.PoissonInterval(0, 2.0, 1e-3, discrete=discrete)


@pytest.mark.parametrize("discrete", [True, False])
def test_binomial_interval(discrete):
    pi = fwdpy11.BinomialInterval(0, 2, 1e-3, discrete=discrete)
    p = pickle.dumps(pi)
    up = pickle.loads(p)
    assert up.beg == pi.beg
    assert up.end == pi.end
    assert up.probability == pi.probability
    assert up.discrete == pi.discrete
    assert up.discrete == discrete

    if discrete is True:
        with pytest.raises(ValueError):
            fwdpy11.BinomialInterval(0, 1, 1e-3, discrete=discrete)
        with pytest.raises(TypeError):
            fwdpy11.BinomialInterval(0.0, 1, 1e-3, discrete=discrete)
        with pytest.raises(TypeError):
            fwdpy11.BinomialInterval(0, 1.0, 1e-3, discrete=discrete)


@pytest.mark.parametrize("discrete", [True, False])
def test_poisson_point(discrete):
    pi = fwdpy11.PoissonPoint(0, 1e-3, discrete=discrete)
    p = pickle.dumps(pi)
    up = pickle.loads(p)
    assert up.position == pi.position
    assert up.mean == pi.mean
    assert up.discrete == pi.discrete
    assert up.discrete == discrete

    if discrete is True:
        with pytest.raises(TypeError):
            fwdpy11.PoissonPoint(0.0, 1e-3, discrete=discrete)


@pytest.mark.parametrize("discrete", [True, False])
def test_binomial_point(discrete):
    pi = fwdpy11.BinomialPoint(0, 1e-3, discrete=discrete)
    p = pickle.dumps(pi)
    up = pickle.loads(p)
    assert up.position == pi.position
    assert up.probability == pi.probability
    assert up.discrete == pi.discrete
    assert up.discrete == discrete

    if discrete is True:
        with pytest.raises(TypeError):
            fwdpy11.BinomialPoint(0.0, 1e-3, discrete=discrete)


@pytest.mark.parametrize("discrete", [True, False])
def test_fixed_crossovers(discrete):
    pi = fwdpy11.FixedCrossovers(0, 2, 3, discrete=discrete)
    p = pickle.dumps(pi)
    up = pickle.loads(p)
    assert up.beg == pi.beg
    assert up.end == pi.end
    assert up.num_xovers == pi.num_xovers
    assert up.discrete == pi.discrete
    assert up.discrete == discrete

    if discrete is True:
        with pytest.raises(ValueError):
            fwdpy11.FixedCrossovers(0, 1, 3, discrete=discrete)
        with pytest.raises(TypeError):
            fwdpy11.FixedCrossovers(0.0, 1, 3, discrete=discrete)
        with pytest.raises(TypeError):
            fwdpy11.FixedCrossovers(0, 1.0, 3, discrete=discrete)
