"""
This module provided pre-calculated demographic models
for human populations.
"""
import enum

import attr
import numpy as np

import fwdpy11

from .demographic_model_details import (DemographicModelCitation,
                                        DemographicModelDetails)


class TennessenModel(enum.Enum):
    """
    Enumeration type declaring alternative implementations
    of the Tennessen et al. (2012) model. See
    :func:`fwdpy11.demographic_models.human.tennessen` for
    details.
    """

    V0 = 0
    """
    Default version of the model
    """
    V1 = 1
    """
    Alternate version of the model
    """


@attr.s()
class _ModelParameterValidator(object):
    """ Make sure that we get an actual int """

    burnin = attr.ib(type=int, kw_only=True, validator=attr.validators.instance_of(int))

    @burnin.validator
    def validate_burnin(self, attribute, value):
        if value < 0:
            raise ValueError("burnin must be >= 0")


def tennessen(burnin: int = 20, model_version: TennessenModel = TennessenModel.V0):
    """
    Generate parameters for the demographic model described in:

    Tennessen, Jacob A., Abigail W. Bigham, Timothy D. O’Connor, Wenqing Fu,
    Eimear E. Kenny, Simon Gravel, Sean McGee, et al. 2012.
    “Evolution and Functional Impact of Rare Coding Variation
    from Deep Sequencing of Human Exomes.” Science 337 (6090): 64–69.

    :param burnin: Burn-in time, as a multiplier of the ancestral population size.
                   Defaults to 20.
    :type param: int
    :param model_version: Specify which "flavor" of the model to run.  See below.
    :type model_version: fwdpy11.demographic_models.human.TennessenModel

    :returns: The demographic model
    :rtype: fwdpy11.demographic_models.DemographicModelDetails

    The ``metadata`` field of the return value contains information about
    the simulation length, ancestral population size (``Nref``), etc.,
    and the mapping of integer values to deme names.

    When this model is run, deme ``0`` corresponds to ``African`` and deme
    ``1`` corresponds to Eurasian.

    Due to uncertainty in the literature over the details of this model,
    we provide two slightly different implementations.  The first
    implementation is the default, and is specified by
    :attr:`fwdpy11.demographic_models.human.TennessenModel.V0`.  For this
    model, the final deme sizes are those shown in Figure 3 of the
    Tennessen et al. (2012) paper.
    If :attr:`fwdpy11.demographic_models.human.TennessenModel.V1` is
    used instead, then the final sizes are obtained by calculating them
    forwards in time using a combination of the information from
    Tennessen et al. (2012) and Fu et al. 2013. However, due to having
    to convert from continuous to discrete time, etc., other definitions
    of this model are also possible.  Here, all deme sizes and times
    are obtained using :func:`numpy.rint`.


    .. note::

        This implementation is based on code provided by Aaron Ragsdale.

    .. versionadded:: 0.8.0
    """
    _ModelParameterValidator(burnin=burnin)
    Nref = 7310  # Ancestral population dize
    NAfr0 = 14474  # Initial size change
    NB = 1861  # Eurasian bottleneck size
    if model_version == TennessenModel.V0:
        NAfr = 420000  # Final African pop'n size
    elif model_version == TennessenModel.V1:
        NAfr = 423125  # Final African pop'n size
    else:
        raise ValueError("invalid model_version")

    NEur0 = 1032  # Eurasian second bottleneck size
    if model_version == TennessenModel.V0:
        NEur1 = 9237  # Eurasian size as rapid growth starts
    elif model_version == TennessenModel.V1:
        NEur1 = 9279  # Eurasian size as rapid growth starts
    else:
        raise ValueError("invalid model_version")

    if model_version == TennessenModel.V0:
        NEur = 512000  # Final Eurasian pop'n size
    elif model_version == TennessenModel.V1:
        NEur = 501425  # Final Eurasian pop'n size
    else:
        raise ValueError("invalid model_version")
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
        f"Sequencing of Human Exomes.” Science 337 (6090): 64–69. ",
        f"and Figure S5 of: ",
        f"Fu, Wenqing, Timothy D. O’Connor, Goo Jun, Hyun Min Kang, Goncalo Abecasis, "
        f"Suzanne M. Leal, Stacey Gabriel, et al. 2013. “Analysis of 6,515 "
        f"Exomes Reveals the Recent Origin of Most Human Protein-Coding Variants.” "
        f"Nature 493 (7431): 216–20",
    )

    full_citation = str().join([str(i) for i in full_citation])

    return DemographicModelDetails(
        model=ddemog,
        name="Tennessen et al. model of African and European demography.",
        source={"function": "fwdpy11.demographic_models.human.tennessen"},
        parameters={"burnin": burnin},
        citation=DemographicModelCitation(
            DOI="10.1126/science.1219240 and 10.1038/nature11690",
            full_citation=full_citation,
            metadata=None,
        ),
        metadata={
            "deme_labels": {0: "African", 1: "Eurasian"},
            "simlen": total_sim_length,
            "Nref": Nref,
        },
    )


def jouganous_three_deme(burnin: int = 20):
    """
    Generate parameters for the demographic model described in Table 2 of:

    Jouganous, Julien, Will Long, Aaron P. Ragsdale, and Simon Gravel. 2017.
    “Inferring the Joint Demographic History of Multiple Populations:
    Beyond the Diffusion Approximation.” Genetics 206 (3): 1549–67.

    :param burnin: Burn-in time, as a multiplier of the ancestral population size.
                   Defaults to 20.
    :type param: int
    :returns: The demographic model
    :rtype: fwdpy11.demographic_models.DemographicModelDetails

    The ``metadata`` field of the return value contains information about
    the simulation length, ancestral population size (``Nref``), etc.,
    and the mapping of integer values to deme names.

    .. versionadded:: 0.8.1

    """
    _ModelParameterValidator(burnin=burnin)

    NA = 11273  # Ancestral population size

    gens_burn_in = burnin * NA

    # TAF years ago, change the ancestral
    # deme size to the African deme size.
    TAF = 312000
    NAF = 23721

    # TB years ago, the ancestral population
    # splits into the population ancestral
    # to modern Africa and population B,
    # which is the ancestor of modern
    # Eurasian populations.
    TB = 125000
    NB = 3104

    # TAS_EU years ago, pop B splits into
    # EU (modern European) and AS (modern
    # Asian).
    TAS_EU = 42300
    NAS0 = 923
    NEU0 = 2271

    generation_time = 29
    # These are the CONTINUOUS growth
    # rates inferred in Jouganous et al.
    rAS = 0.00309  # units of 2NA
    rEU = 0.00196
    tg = TAS_EU / generation_time  # Time AS and EU grow exponentially, in generations

    # Get the discrete-time growth rates...
    NASF = NAS0 * np.exp(rAS * tg)
    NEUF = NEU0 * np.exp(rEU * tg)
    NASFi = np.rint(NASF)
    NEUFi = np.rint(NEUF)

    # ... which are these:
    dGAS = np.exp((np.log(NASFi) - np.log(NAS0)) / np.rint(tg))
    dGEU = np.exp((np.log(NEUFi) - np.log(NEU0)) / np.rint(tg))

    # Now, convert all the times from years into discrete
    # generations.

    TAFg = 0

    TBg = np.rint((TAF - TB) / generation_time).astype(int)
    TAS_EUg = np.rint((TB - TAS_EU) / generation_time).astype(int)
    total_sim_length = gens_burn_in + TAFg + TBg + TAS_EUg + np.rint(tg).astype(int)

    # Build up the events list
    deme_labels = {"AFR": 0, "B_CEU": 1, "CHB": 2}

    migmatrix = np.zeros(9).reshape(3, 3)
    migmatrix[0, 0] = 1.0
    change_migrates = []

    T = TAFg

    # Change size of AFR
    set_deme_sizes = [
        fwdpy11.SetDemeSize(
            when=gens_burn_in + T, deme=deme_labels["AFR"], new_size=NAF
        )
    ]

    T += TBg

    # "Bud" off deme B
    mass_migrations = [
        fwdpy11.copy_individuals(
            when=gens_burn_in + T,
            source=deme_labels["AFR"],
            destination=deme_labels["B_CEU"],
            fraction=NB / NAF,
        )
    ]

    mAF_B = 15.8e-5
    change_migrates.append(
        fwdpy11.SetMigrationRates(
            when=gens_burn_in + T,
            deme=deme_labels["AFR"],
            migrates=[1.0 - mAF_B, mAF_B, 0.0],
        )
    )
    change_migrates.append(
        fwdpy11.SetMigrationRates(
            when=gens_burn_in + T,
            deme=deme_labels["B_CEU"],
            migrates=[mAF_B, 1.0 - mAF_B, 0.0],
        )
    )

    # Split B into CEU and CHB ("European" and "Asian")
    # CEU inherits the deme label of B.
    T += TAS_EUg
    mass_migrations.append(
        fwdpy11.copy_individuals(
            when=gens_burn_in + T,
            source=deme_labels["B_CEU"],
            destination=deme_labels["CHB"],
            fraction=NAS0 / NB,
        )
    )

    set_deme_sizes.append(
        fwdpy11.SetDemeSize(
            when=gens_burn_in + T, deme=deme_labels["B_CEU"], new_size=NEU0
        )
    )

    set_growth_rates = [
        fwdpy11.SetExponentialGrowth(
            when=gens_burn_in + T, deme=deme_labels["B_CEU"], G=dGEU
        )
    ]

    set_growth_rates.append(
        fwdpy11.SetExponentialGrowth(
            when=gens_burn_in + T, deme=deme_labels["CHB"], G=dGAS
        )
    )

    mAF_CEU = 1.10e-5
    mAF_CHB = 0.48e-5
    mCEU_CHB = 4.19e-5

    change_migrates.append(
        fwdpy11.SetMigrationRates(
            when=gens_burn_in + T,
            deme=deme_labels["AFR"],
            migrates=[1.0 - (mAF_CEU + mAF_CHB), mAF_CEU, mAF_CHB],
        )
    )

    change_migrates.append(
        fwdpy11.SetMigrationRates(
            when=gens_burn_in + T,
            deme=deme_labels["B_CEU"],
            migrates=[mAF_CEU, 1.0 - (mAF_CEU + mCEU_CHB), mCEU_CHB],
        )
    )

    change_migrates.append(
        fwdpy11.SetMigrationRates(
            when=gens_burn_in + T,
            deme=deme_labels["CHB"],
            migrates=[mAF_CHB, mCEU_CHB, 1.0 - (mAF_CHB + mCEU_CHB)],
        )
    )

    ddemog = fwdpy11.DiscreteDemography(
        mass_migrations=mass_migrations,
        set_deme_sizes=set_deme_sizes,
        migmatrix=migmatrix,
        set_migration_rates=change_migrates,
        set_growth_rates=set_growth_rates,
    )

    full_citation = (
        f"Jouganous, Julien, Will Long, Aaron P. Ragsdale, and Simon Gravel. 2017."
        f" “Inferring the Joint Demographic History "
        f"of Multiple Populations: Beyond the Diffusion Approximation.”"
        f" Genetics 206 (3): 1549–67."
    )
    full_citation = str().join([str(i) for i in full_citation])

    return DemographicModelDetails(
        model=ddemog,
        name="Jouganous et al. 3-deme model.",
        source={"function": "fwdpy11.demographic_models.human.jouganous_three_deme"},
        parameters={"burnin": burnin},
        citation=DemographicModelCitation(
            DOI="10.1534/genetics.117.200493",
            full_citation=full_citation,
            metadata="Parameters from Table 2.",
        ),
        metadata={
            "deme_labels": {value: key for key, value in deme_labels.items()},
            "simlen": total_sim_length,
            "Nref": NA,
        },
    )
