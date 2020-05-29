"""
This module provided pre-calculated demographic models
for human populations.
"""
import attr
import numpy as np

import fwdpy11

from .demographic_model_details import (DemographicModelCitation,
                                        DemographicModelDetails)


@attr.s()
class _TennessenParameterValidator(object):
    """ Make sure that we get an actual int """

    burnin = attr.ib(type=int, validator=attr.validators.instance_of(int))


def tennessen(burnin: int = 20):
    """
    Generate parameters for the demographic model described in:

    Tennessen, Jacob A., Abigail W. Bigham, Timothy D. O’Connor, Wenqing Fu,
    Eimear E. Kenny, Simon Gravel, Sean McGee, et al. 2012.
    “Evolution and Functional Impact of Rare Coding Variation
    from Deep Sequencing of Human Exomes.” Science 337 (6090): 64–69.

    :param burnin: Burn-in time, as a multiplier of the ancestral population size.
                   Defaults to 20.
    :type param: int

    :returns: The demographic model
    :rtype: fwdpy11.demographic_models.DemographicModelDetails

    The ``metadata`` field of the return value contains information about
    the simulation length, ancestral population size (``Nref``), etc.,
    and the mapping of integer values to deme names.

    When this model is run, deme ``0`` corresponds to ``African`` and deme
    ``1`` corresponds to Eurasian.

    .. note::

        This implementation is based on code provided by Aaron Ragsdale.

    .. versionadded:: 0.7.2

    .. versionchanged:: 0.8.0

        Returns instance of :class:`fwdpy11.demographic_models.DemographicModelDetails`

    """
    _TennessenParameterValidator(burnin)
    Nref = 7310  # Ancestral population dize
    NAfr0 = 14474  # Initial size change
    NB = 1861  # Eurasian bottleneck size
    NAfr = 424000  # Final African pop'n size
    NEur0 = 1032  # Eurasian second bottleneck size
    NEur1 = 9237  # Eurasian size as rapid growth starts
    NEur = 512000  # Final Eurasian pop'n size
    T_Af = 148000
    T_B = 51000
    T_Eu_As = 23000
    T_accel = 5115
    mB = 15e-5
    mF = 2.5e-5
    generation_time = 25

    # List of demographic events:
    # keep track of size change, copying, and migration rate change events in
    # separate lists
    size_change = []
    copy = []
    mig_rates = []
    growth_rates = []

    # number of generations in epochs
    T0 = np.rint((T_Af - T_B) / generation_time).astype(int)  # pre-split
    # split to bottleneck, no growth
    T1 = np.rint((T_B - T_Eu_As) / generation_time).astype(int)
    T2 = np.rint((T_Eu_As - T_accel) / generation_time).astype(
        int
    )  # Eu growth with r_Eu0
    # accelerated growth in Af and Eu
    T3 = np.rint(T_accel / generation_time).astype(int)

    M_init = np.zeros(4).reshape(2, 2)
    M_init[0, 0] = 1
    mm = fwdpy11.MigrationMatrix(M_init)

    # burn in for 20*Ne generations
    gens_burn_in = burnin * Nref
    total_sim_length = gens_burn_in + T0 + T1 + T2 + T3

    # init: size change of common ancestral population
    size_change.append(fwdpy11.SetDemeSize(when=gens_burn_in, deme=0, new_size=NAfr0))

    # T0: mass migration, copy from A to Eu bottleneck population
    copy.append(
        fwdpy11.copy_individuals(
            when=gens_burn_in + T0, source=0, destination=1, fraction=NB / NAfr0,
        )
    )
    size_change.append(fwdpy11.SetDemeSize(when=gens_burn_in + T0, deme=1, new_size=NB))
    # at the same time, set migration rate between deme 0 and 1 to m_A_B
    mig_rates.append(fwdpy11.SetMigrationRates(gens_burn_in + T0, 0, [1 - mB, mB]))
    mig_rates.append(fwdpy11.SetMigrationRates(gens_burn_in + T0, 1, [mB, 1 - mB]))

    # T1: adjust size of Eu to Eu0 and set growth rate
    size_change.append(
        fwdpy11.SetDemeSize(when=gens_burn_in + T0 + T1, deme=1, new_size=NEur0)
    )
    r_Eur0 = (NEur1 / NEur0) ** (1 / T2) - 1
    growth_rates.append(
        fwdpy11.SetExponentialGrowth(when=gens_burn_in + T0 + T1, deme=1, G=1 + r_Eur0)
    )
    # set migration rates to contemporary rates
    mig_rates.append(fwdpy11.SetMigrationRates(gens_burn_in + T0 + T1, 0, [1 - mF, mF]))
    mig_rates.append(fwdpy11.SetMigrationRates(gens_burn_in + T0 + T1, 1, [mF, 1 - mF]))

    # T2: set growth rates to accelerated rates in both populations
    r_AfrF = (NAfr / NAfr0) ** (1 / T3) - 1
    r_EurF = (NEur / NEur1) ** (1 / T3) - 1
    growth_rates.append(
        fwdpy11.SetExponentialGrowth(
            when=gens_burn_in + T0 + T1 + T2, deme=0, G=1 + r_AfrF
        )
    )
    growth_rates.append(
        fwdpy11.SetExponentialGrowth(
            when=gens_burn_in + T0 + T1 + T2, deme=1, G=1 + r_EurF
        )
    )

    ddemog = fwdpy11.DiscreteDemography(
        mass_migrations=copy,
        set_deme_sizes=size_change,
        migmatrix=mm,
        set_migration_rates=mig_rates,
        set_growth_rates=growth_rates,
    )
    full_citation = (
        f"Tennessen, Jacob A., Abigail W. Bigham, Timothy D. O’Connor, "
        f"Wenqing Fu, Eimear E. Kenny, Simon Gravel, Sean McGee, et al. 2012. "
        f"“Evolution and Functional Impact of Rare Coding Variation from Deep "
        f"Sequencing of Human Exomes.” Science 337 (6090): 64–69.",
    )
    return DemographicModelDetails(
        model=ddemog,
        name="Tennessen et al. model of African and European demography.",
        source={"function": "fwdpy11.demographic_models.human.tennessen"},
        parameters={"burnin": burnin},
        citation=DemographicModelCitation(
            DOI="10.1126/science.1219240", full_citation=full_citation[0], metadata=None
        ),
        metadata={
            "deme_labels": {0: "African", 1: "Eurasian"},
            "simlen": total_sim_length,
            "Nref": Nref,
        },
    )
