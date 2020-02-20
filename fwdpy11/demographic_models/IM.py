"""
This module provides functions to generate demographic events for
"isolation-with-migration", or IM, models.
"""


def two_deme_IM(Nanc, T, psplit, Ns, migrates, burnin=10.0,
                as_dict=False):
    """
    Isolation-with-migration (IM) model for two demes.

    An ancestral population splits into two daughter demes.
    At the time of the split, ``psplit`` of the ancestral
    population moves into deme 1.  The two daughter populations
    begin exponential growth until the present time and migration
    may occur between them.

    :param Nanc: The ancestral population size.
    :type Nanc: int
    :param T: The time of the split, in units of Nanc generations
              into the past.
    :type T: float
    :param psplit: The proportion of the ancestral population that splits
                   off to found deme 1
    :type psplit: float
    :param Ns: The final sizes of demes 0 and 1, relative to Nanc
    :type Ns: tuple
    :param migrates: The migration rates from 0 to 1 and from 1 to 0,
                     respectively. Migration rates are the fraction
                     of the destination deme replaced by the source
                     deme.
    :type migrates: float
    :param burnin: Time to simulate before the split, in units of Nanc
    :type burnin: float
    :param as_dict: If False, return an instance of
                    :class:`fwdpy11.DiscreteDemography`.
                    Otherwise, return a dict that can be
                    used to construct a ``DiscreteDemography``
                    instance.
    :type as_dict: boolean
    :rtype: tuple
    :returns: the model events, the split time, and
              the time to simulate post-split.

    .. note::

        The events returned by this model assume/require that you will
        construct a population with intitial size ``Nanc``.

    .. versionadded:: 0.6.0
    """
    import fwdpy11
    import numpy as np
    N0, N1 = Ns
    m01, m10 = migrates

    split_time = np.rint(Nanc*burnin).astype(int)
    # The split event
    split = [fwdpy11.move_individuals(when=split_time, source=0,
                                      destination=1,
                                      fraction=psplit)]
    # Get growth rates and set growth rate changes,
    # taking care to handle our rounding!
    gens_post_split = np.rint(Nanc*T).astype(int)
    N0split = np.rint(Nanc*(1.-psplit))
    if N0split == 0 or N0split == Nanc:
        raise ValueError("invalid value for psplit: {}".format(psplit))
    N0final = np.rint(N0*Nanc)
    N1split = np.rint(Nanc*psplit)
    if N1split == 0 or N1split == Nanc:
        raise ValueError("invalid value for psplit: {}".format(psplit))
    N1final = np.rint(N1*Nanc)
    G0 = fwdpy11.exponential_growth_rate(N0split, N0final,
                                         gens_post_split)
    G1 = fwdpy11.exponential_growth_rate(N1split, N1final,
                                         gens_post_split)
    growth = [fwdpy11.SetExponentialGrowth(split_time, 0, G0),
              fwdpy11.SetExponentialGrowth(split_time, 1, G1)]

    # Set up the migration matrix for two demes, but only
    # deme zero exists.
    m = fwdpy11.migration_matrix_single_extant_deme(2, 0)
    # The rows of the matrix change at the split:
    cm = [fwdpy11.SetMigrationRates(split_time, 0, [1.-m10, m10]),
          fwdpy11.SetMigrationRates(split_time, 1, [m01, 1.-m01])]

    mdict = {'mass_migrations': split,
             'set_growth_rates': growth,
             'set_migration_rates': cm,
             'migmatrix': m
             }
    if as_dict is False:
        return fwdpy11.DiscreteDemography(**mdict), split_time, gens_post_split

    return mdict, split_time, gens_post_split
