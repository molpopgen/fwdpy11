from .wfevolve import evolve_singlepop_regions_cpp

def evolve_regions_sampler_fitness(rng,pop,popsizes,mu_neutral,
        mu_selected,recrate,nregions,sregions,recregions,fitness,
        recorder,selfing_rate = 0.):
    from ..regions import makeMutationRegions,makeRecombinationRegions
    mm=makeMutationRegions(nregions,sregions)
    rm=makeRecombinationRegions(recregions)
    evolve_singlepop_regions_cpp(rng,pop,popsizes,mu_neutral,
            mu_selected,recrate,mm,rm,fitness,recorder,selfing_rate)
