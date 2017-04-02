from .wfevolve_qtrait import evolve_singlepop_regions_qtrait_cpp
def evolve_regions_sampler_fitness(rng,pop,popsizes,mu_neutral,
        mu_selected,recrate,nregions,sregions,recregions,trait_model,
        recorder,environments,selfing_rate = 0.):
    """
    Evolve a single deme according to a Wright-Fisher life cycle 
    with arbitrary changes in population size, a specified fitness model,
    and a temporal sampler.

    The fitness model is Gaussian stabilizing selection with respect to an
    optimum:

    .. math::

        w = e^{-\\frac{(P-O)^2}{2VS}}

    The trait value is :math:`P=G+E`, where :math:`G` is the genetic component, calculated
    using an object of type `fwdpy11.fitness.SpopFitness`.  The "environmental" term, 
    :math:`E` is a Gassian deviate with mean 0 and standard deviation 
    :math:`\sigma_e`.

    :param rng: A :class:`fwdpy11.fwdpy11_types.GSLrng`
    :param pop: A :class:`fwdpy11.fwdpy11_types.Spop`
    :param popsizes: A 1d NumPy array representing population sizes over time.
    :param mu_neutral: The neutral mutation rate (per gamete, per generation)
    :param mu_selected: The selected mutation rate (per gamete, per generation)
    :param recrate: The recombination reate (per diploid, per generation)
    :param nregions: A list of :class:`fwdpy11.regions.Region`.
    :param sregions: A list of :class:`fwdpy11.regions.Sregion`.
    :param recregions: A list of :class:`fwdpy11.regions.Region`.
    :param trait_model: A :class:`fwdpy11.fitness.SpopFitness`.
    :param recorder: A callable to record data from the population.
    :param environments: A list of environmental conditions (see below).
    :param selfing_rate: (default 0.0) The probability than an individual selfs.

    The environmental conditions are specified by a list of tuples containing
    four elements each.  These elements are:

    * Generation.  
    * optimum, the optimum trait value.
    * VS, the intensity of selection against deviations from the optimum trait value.
    * :math:`\sigma_e`, the standard deviation of random noise added to trait value.

    When the population's generation first becomes :math:`\geq` `Generation`, then 
    the values for `optimum`, `VS` and :math:`\sigma_e` are applied.
    """
    from .internal import makeMutationRegions,makeRecombinationRegions
    mm=makeMutationRegions(nregions,sregions)
    rm=makeRecombinationRegions(recregions)
    evolve_singlepop_regions_qtrait_cpp(rng,pop,popsizes,mu_neutral,
            mu_selected,recrate,mm,rm,trait_model,recorder,selfing_rate,environments)
