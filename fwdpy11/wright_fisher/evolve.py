from .wfevolve import evolve_singlepop_regions_cpp
#void
#evolve_singlepop_regions_cpp(
#    const fwdpy11::GSLrng_t& rng, fwdpy11::singlepop_t& pop, const unsigned& N,
#    const unsigned generations, const double mu_neutral,
#    const double mu_selected, const double recrate,
#    const KTfwd::extensions::discrete_mut_model& mmodel,
#    const KTfwd::extensions::discrete_rec_model& rmodel,
#    const fwdpy11::singlepop_fitness& fitness,
#    fwdpy11::singlepop_temporal_sampler recorder, const double selfing_rate)

def evolve_regions_sampler_fitness(rng,pop,N,generations,mu_neutral,
        mu_selected,recrate,nregions,sregions,recregions,fitness,
        recorder,selfing_rate = 0.):
    from ..regions import makeMutationRegions,makeRecombinationRegions
    mm=makeMutationRegions(nregions,sregions)
    rm=makeRecombinationRegions(recregions)
    evolve_singlepop_regions_cpp(rng,pop,N,generations,mu_neutral,
            mu_selected,recrate,mm,rm,fitness,recorder,selfing_rate)
