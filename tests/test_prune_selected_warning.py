import pytest

import fwdpy11


@pytest.mark.parametrize(
    "gvalue",
    [
        fwdpy11.Additive(scaling=2.0),
        fwdpy11.GBR(
            fwdpy11.GaussianStabilizingSelection.single_trait(
                [fwdpy11.Optimum(when=0, optimum=1.0, VS=1.0)]
            )
        ),
        fwdpy11.Multiplicative(
            scaling=2,
            gvalue_to_fitness=fwdpy11.GaussianStabilizingSelection.single_trait(
                [fwdpy11.Optimum(when=0, optimum=1.0, VS=1.0)]
            ),
        ),
        fwdpy11.Multiplicative(
            scaling=2,
            noise=fwdpy11.GaussianNoise(mean=0.0, sd=1.0),
        ),
    ],
)
def test_warning(gvalue):
    pdict = {
        "gvalue": gvalue,
        "rates": (0, 0, 0),
        "demography": fwdpy11.ForwardDemesGraph.tubes([100], burnin=1),
        "simlen": 100,
    }
    with pytest.warns(fwdpy11.PruneSelectedWarning):
        params = fwdpy11.ModelParams(**pdict)
        assert params.prune_selected is True


@pytest.mark.parametrize(
    "gvalue",
    [
        fwdpy11.Additive(scaling=2.0),
        fwdpy11.GBR(
            fwdpy11.GaussianStabilizingSelection.single_trait(
                [fwdpy11.Optimum(when=0, optimum=1.0, VS=1.0)]
            )
        ),
        fwdpy11.Multiplicative(
            scaling=2,
            gvalue_to_fitness=fwdpy11.GaussianStabilizingSelection.single_trait(
                [fwdpy11.Optimum(when=0, optimum=1.0, VS=1.0)]
            ),
        ),
        fwdpy11.Multiplicative(
            scaling=2,
            noise=fwdpy11.GaussianNoise(mean=0.0, sd=1.0),
        ),
    ],
)
def test_no_warning(gvalue):
    pdict = {
        "gvalue": gvalue,
        "rates": (0, 0, 0),
        "demography": fwdpy11.ForwardDemesGraph.tubes([100], burnin=1),
        "simlen": 100,
        "prune_selected": False,
    }

    try:
        params = fwdpy11.ModelParams(**pdict)
        assert params.prune_selected is False
    except fwdpy11.PruneSelectedWarning as e:
        pytest.fail(f"{e}")


@pytest.mark.parametrize(
    "gvalue",
    [
        fwdpy11.Multiplicative(
            scaling=2,
            noise=fwdpy11.GaussianNoise(mean=0.0, sd=1.0),
        ),
    ],
)
@pytest.mark.parametrize("prune_selected", [False, True])
def test_multiplicative_fitness_model(gvalue, prune_selected):
    pdict = {
        "gvalue": fwdpy11.Multiplicative(scaling=2.0),
        "rates": (0, 0, 0),
        "demography": fwdpy11.ForwardDemesGraph.tubes([100], burnin=1),
        "simlen": 100,
        "prune_selected": prune_selected,
    }
    try:
        params = fwdpy11.ModelParams(**pdict)
        assert params.prune_selected == prune_selected
    except fwdpy11.PruneSelectedWarning as e:
        pytest.fail(f"unexpected exception {e}")
