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

There are a number of equivalent ways to do this.  They differ in efficiency.

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

The most straightforward way to do this is with collections.Counter.

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
zip:

.. testcode::

    #Filter out extinct mutations and neutral mutations in the list comprehension
    freq_esize = np.array([(float(i)/float(2*pop.N),j.s) for i,j in zip(pop.mcounts,pop.mutations) 
        if i > 0 and j.neutral is False],dtype=[('freq',np.float64),('esize',np.float64)])


    

