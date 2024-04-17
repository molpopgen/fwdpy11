import numpy as np
import pytest

import fwdpy11


def test_single_deme_moving_optimum():
    # Copied from test_tree_sequences.pyh::set_up_quant_trait_model,
    # but the first "when" is edited to trigger an error
    N = 1000
    demography = fwdpy11.ForwardDemesGraph.tubes(
        [N], int(np.rint(N)), burnin_is_exact=True
    )
    rho = 2.0
    r = rho / (4 * N)
    Opt = fwdpy11.Optimum
    GSSmo = fwdpy11.GaussianStabilizingSelection.single_trait(
        [Opt(when=10, optimum=0.0, VS=1.0), Opt(when=N, optimum=1.0, VS=1.0)]
    )
    a = fwdpy11.Additive(2.0, GSSmo)
    p = {
        "nregions": [],
        "sregions": [fwdpy11.GaussianS(0, 1, 1, 0.25)],
        "recregions": [fwdpy11.PoissonInterval(0, 1, r)],
        "rates": (0.0, 0.025, None),
        "gvalue": a,
        "prune_selected": False,
        "demography": demography,
        "simlen": np.rint(N).astype(int),
    }
    with pytest.raises(fwdpy11.TimingError):
        _ = fwdpy11.ModelParams(**p)


def test_two_deme_moving_optimum():
    # Example 07 from the demes tutorial
    yaml = """
    time_units: generations
    demes:
      - name: X
        epochs:
          - end_time: 1000
            start_size: 2000
      - name: A
        ancestors:
          - X
        epochs:
          - start_size: 2000
      - name: B
        ancestors:
          - X
        epochs:
          - start_size: 2000
    """
    demography = fwdpy11.ForwardDemesGraph.from_demes(
        yaml,
        burnin=100,
        burnin_is_exact=True,
    )
    Opt = fwdpy11.Optimum
    GSSmo_ancestor = fwdpy11.GaussianStabilizingSelection.single_trait(
        [Opt(when=0, optimum=0.0, VS=1.0), Opt(when=100, optimum=1.0, VS=1.0)]
    )
    gvalue_ancestor = fwdpy11.Additive(2.0, GSSmo_ancestor)
    GSSmo_daughter_1 = fwdpy11.GaussianStabilizingSelection.single_trait(
        [Opt(when=100, optimum=1.0, VS=1.0)]
    )
    gvalue_daughter_1 = fwdpy11.Additive(2.0, GSSmo_daughter_1)

    # This is the error: the first "when" value is set AFTER
    # this population first exists, resulting in a "late"
    # catch at run time in evolvets
    GSSmo_daughter_2 = fwdpy11.GaussianStabilizingSelection.single_trait(
        [Opt(when=110, optimum=1.0, VS=1.0)]
    )
    gvalue_daughter_2 = fwdpy11.Additive(2.0, GSSmo_daughter_2)
    p = {
        "nregions": [],
        "sregions": [fwdpy11.GaussianS(0, 1, 1, 0.25)],
        "recregions": [fwdpy11.PoissonInterval(0, 1, 1e-3)],
        "rates": (0.0, 0.025, None),
        "gvalue": [gvalue_ancestor, gvalue_daughter_1, gvalue_daughter_2],
        "prune_selected": False,
        "demography": demography,
        "simlen": demography.final_generation,
    }
    with pytest.raises(fwdpy11.TimingError):
        _ = fwdpy11.ModelParams(**p)
