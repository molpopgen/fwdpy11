.. _processingpops:

Processing simulated populations
======================================================================

We will work through some basic operations on populations.  The main point is to show how one can access the data
structures described in :ref:`data_types`.  

You should read the following sections first:

* :ref:`data_types`
* :ref:`model_params`

First, we'll quickly simulate a single deme for `N` generations:

.. testcode::

    import fwdpy11
    import fwdpy11.fitness
    import fwdpy11.model_params
    import fwdpy11.ezparams
    import fwdpy11.wright_fisher
    import numpy as np

    rng = fwdpy11.GSLrng(42)

    theta,rho = 100.0,100.0
    pop = fwdpy11.SlocusPop(1000)

    pdict = fwdpy11.ezparams.mslike(pop,simlen=pop.N,
        dfe=fwdpy11.regions.ExpS(0,1,1,-0.1,1),pneutral = 0.95)
    params = fwdpy11.model_params.SlocusParams(**pdict)

    fwdpy11.wright_fisher.evolve(rng,pop,params)

Population mean fitness
------------------------------------------------------------------------------------------------------

There are a number of equivalent ways to do this.  They differ in efficiency.  The fitness is stored
as the `w` field of a diploid.

Method 1 is a direct summation:

.. testcode::

    wbar1 = 0.0
    for i in pop.diploids:
        wbar1 += i.w
    wbar1 /= float(pop.N)

Method 2 involves a sum over a list comprehension over diploids:

.. testcode::

    wbar2 = sum([i.w for i in pop.diploids])/float(pop.N)

Method 3 uses Numpy because we can:

.. testcode::

    import numpy as np
    wbar3 = np.array([i.w for i in pop.diploids]).sum()/float(pop.N)

The three methods are all conceptually identical and will be numerically identical modulo rounding 
errors due to non-transitivity of floating-point addition.

The site-frequency spectrum at neutral and non-neutral positions
------------------------------------------------------------------------------------------------------

The most straightforward way to do this is with collections.Counter. We need to iterate over `mcounts`
to filter out extinct mutations and also test that each *extant* variant is non-neutral.

.. testcode::

    import collections
    sfsn = collections.Counter()   
    sfss = collections.Counter()   

    for i in range(len(pop.mcounts)):
        #Skip extinct variants:
        if pop.mcounts[i] > 0 :  
            #Distinguish neutral from non-neutral
            #mutations
            if pop.mutations[i].neutral is False:
                sfss[pop.mcounts[i]] += 1
            else:
                sfsn[pop.mcounts[i]] += 1

The relationship between frequency and effect size
------------------------------------------------------------------------------------------------------

Let's store the result in a numpy structured array.  We can do it in a one-liner involving a list comprehension and a
zip over the `mcounts` and `mutations`:

.. testcode::

    #Filter out extinct mutations and neutral mutations in the list comprehension
    freq_esize = np.array([(float(i)/float(2*pop.N),j.s) for i,j in zip(pop.mcounts,pop.mutations) 
        if i > 0 and j.neutral is False],dtype=[('freq',np.float64),('esize',np.float64)])

Mean number of selected mutations per diploid
------------------------------------------------------------------------------------------------------

The `first` and `second` properties of a diploid are the indexes to that diploid's gametes.  We need
those indexes to access the length of the `smutations` property of each gamete.

.. testcode::

    nselected_per_dip = np.array([len(pop.gametes[i.first].smutations) + len(pop.gametes[i.second].smutations)
        for i in pop.diploids])
    mean_selected_muts_per_diploid = nselected_per_dip.mean()
    
Sum of effect sizes on each gamete in each diploid
------------------------------------------------------------------------------------------------------

We will break the calculation into two steps:

1. Calculate the sum of effect sizes on each extant gamete and record the result in a list mapping gamete index to the
   sum.
2. Iterate over the diploids and use the results from step 1 to get our answer

.. testcode::
    
    def sum_esizes_gamete(pop,i):
        return sum([pop.mutations[m].s for m in pop.gametes[i].smutations])
    #Step 1
    gamete_sum_esizes = {i:sum_esizes_gamete(pop,i) for i in range(len(pop.gametes)) if pop.gametes[i].n > 0}
    #Step 2
    sum_ezizes_per_dip =[(gamete_sum_esizes[i.first],gamete_sum_esizes[i.second]) for i in pop.diploids]

Get a random set of diploids
------------------------------------------------------------------------------------------------------

A common analysis involves a random set of diploids.  For example, a sample of size :math:`n \ll N` is typically assumed
in much of standard coalescent theory.  One could take a random contiguous slice by sampling a random initial offset,
taking care to avoid going past the end of the container (remember that the opaque lists do not support negative
indexing, etc.!).  A perhaps less error-prone method involves just asking `numpy` to do the work for us:

.. testcode::

    #Number of diploids to sample
    n = 50
    np.random.seed(42)
    dips = np.random.choice(pop.N,n,replace=False)
    sample = [pop.diploids[i] for i in dips]
    #Our new container is a Python list.
    #If we had used slices, it would be
    #an opaque list
    assert(type(sample) is list)
    assert(type(sample[0]) is fwdpy11.fwdpy11_types.SingleLocusDiploid)



