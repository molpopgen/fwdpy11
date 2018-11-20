.. _advancedmanip:

Various ways to manipulate a population
======================================================================

This is an example of various ways to extract data from a population.  It complements :ref:`processingpops` and also 
gives me a place to put examples based on specific questions from users.

.. ipython:: python

    import fwdpy11 as fp11
    import fwdpy11.genetic_values
    import fwdpy11.wright_fisher as fp11wf
    import fwdpy11.model_params
    import numpy as np
    import gzip,bz2,lzma
    import fwdpy11.sampling,pickle,lzma,gzip
    N, theta, rho = 1000, 10,10 

    nlist=np.array([N]*(N),dtype=np.uint32)
    nregions = [fp11.Region(beg=1,end=10000,weight=1)]
    sregions=[fp11.ExpS(1,10000,1,-0.25)]
    recregions = [fp11.Region(beg=1,end=10000,weight=1)]
    rng = fp11.GSLrng(42)

    pdict = {'nregions':nregions,
            'sregions':sregions,
            'recregions':recregions,
            'demography':nlist,
            'rates':(theta/float(4*N),1e-2,rho/float(4*N)),
            'gvalue':fwdpy11.genetic_values.SlocusMult(2.)
            }
    params = fwdpy11.model_params.ModelParams(**pdict)

    pop = fp11.SlocusPop(N)
    pops = fp11wf.evolve(rng, pop, params)

Exracting info about the mutations
----------------------------------------------------------------

Get a tuple of the mutation objects
and their counts in the population.
We have to filter out extinct mutations here,
which fwdpp keeps around for memory recycling.

.. ipython:: python

    extant_muts = [(i,j) for i,j in zip(pop.mutations,pop.mcounts) if j > 0]

Print mutation position,s,h,origin time:

.. ipython:: python

    print(str(extant_muts[0][0]))

Exracting info about the diploids
----------------------------------------------------------------

Let's record the position, effect size, and origin time for each mutation in each diploid in a 
NumPy structured array.  We'll do it just for the first gamete in each diploid:

.. ipython:: python

    for dip in pop.diploids:
        gam1_data = np.array([(pop.mutations[i].pos,
            pop.mutations[i].s,
            pop.mutations[i].g) for i in pop.gametes[dip.first].smutations],
            dtype=[('pos',np.float),('s',np.float),('g',np.uint32)]) 

Population mean fitness:

.. ipython:: python

    w = np.array([i.w for i in pop.diploid_metadata])
    print(w.mean())

fwdpy11's fitness/trait value calculators can be called from outside of a running simulation:

.. ipython:: python

    a=fwdpy11.genetic_values.SlocusMult(2.0)

    print(pop.diploid_metadata[0].w,a(0, pop))
    for i in [(pop.diploid_metadata[i].w,a(i,pop)) for i in range(pop.N)]:
        assert(i[0] == i[1])
