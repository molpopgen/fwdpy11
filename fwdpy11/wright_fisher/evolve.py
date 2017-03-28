from .wfevolve import evolve_singlepop_regions_cpp

def evolve_regions_sampler_fitness(rng,pop,popsizes,mu_neutral,
        mu_selected,recrate,nregions,sregions,recregions,fitness,
        recorder,selfing_rate = 0.):
    """
    Evolve a single deme according to a Wright-Fisher life cycle 
    with arbitrary changes in population size, a specified fitness model,
    and a temporal sampler.

    :param rng: A :class:`fwdpy11.fwdpy11_types.GSLrng`
    :param pop: A :class:`fwdpy11.fwdpy11_types.Spop`
    :param popsizes: A 1d NumPy array representing population sizes over time.
    :param mu_neutral: The neutral mutation rate (per gamete, per generation)
    :param mu_selected: The selected mutation rate (per gamete, per generation)
    :param recrate: The recombination reate (per diploid, per generation)
    :param nregions: A list of :class:`fwdpy11.regions.Region`.
    :param sregions: A list of :class:`fwdpy11.regions.Sregion`.
    :param recregions: A list of :class:`fwdpy11.regions.Region`.
    :param fitness: A :class:`fwdpy11.fwdpy11.
    :param recorder: A callable to record data from the population.
    :param selfing_rate: (default 0.0) The probability than an individual selfs.
    """
    from ..internal import makeMutationRegions,makeRecombinationRegions
    mm=makeMutationRegions(nregions,sregions)
    rm=makeRecombinationRegions(recregions)
    evolve_singlepop_regions_cpp(rng,pop,popsizes,mu_neutral,
            mu_selected,recrate,mm,rm,fitness,recorder,selfing_rate)
