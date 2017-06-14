.. _processingpopsNP:

Processing simulated populations using fwdpy11's NumPy interface
======================================================================

You should read the following sections first:

* :ref:`processingpops`

First, we'll quickly simulate a single deme for `N` generations:

.. ipython:: python

    import fwdpy11
    import fwdpy11.fitness
    import fwdpy11.model_params
    import fwdpy11.ezparams
    import fwdpy11.wright_fisher
    import numpy as np

    rng = fwdpy11.GSLrng(42)

    theta,rho = 100.0,100.0
    pop = fwdpy11.SlocusPop(1000)

    pdict = fwdpy11.ezparams.mslike(pop,simlen=pop.N,dfe=fwdpy11.regions.ExpS(0,1,1,-0.1,1),pneutral = 0.95)

    params = fwdpy11.model_params.SlocusParams(**pdict)
    fwdpy11.wright_fisher.evolve(rng,pop,params)

Accessing diploid genetic values and fitness
--------------------------------------------------------------------------------------

To get the g,e, and w fields of a diploid:

.. ipython:: python

    dip_traits = np.array(pop.diploids.trait_array())
    dip_traits.dtype

You can easily get mean fitness, etc., now:

.. ipython:: python

    #Use our structured array:
    dip_traits['w'].mean()
    #Access via the diploid objects + list comprehension:
    np.array([i.w for i in pop.diploids]).mean()
    #Do it with pure Python:
    sum([i.w for i in pop.diploids])/float(len(pop.diploids))

.. note::
    The above three calculations give two different numerical results.
    Both are, however, correct given the limits of floating point arithmetic.

Performance
++++++++++++++++++++++++++++++++++++++++++++++++++++++

What's the fastest way to get mean fitness?

.. ipython:: python

    %timeit np.array(pop.diploids.trait_array())['w'].mean()

.. ipython:: python

    %timeit np.array([i.w for i in pop.diploids]).mean()

.. ipython:: python

    %timeit sum([i.w for i in pop.diploids])/float(len(pop.diploids))

Using our stuctured array is vastly more efficient.  For the record, a simple Python loop is the slowest way to go:

.. ipython:: python

    def wbar(pop):
        s=0.0
        for i in range(len(pop.diploids)):
            s += pop.diploids[i].w
        return s/float(len(pop.diploids))

    %timeit wbar(pop)


Slicing, etc.
++++++++++++++++++++++++++++++++++++++++++++++++++++++

To get these fields for a subset of individuals, index with a numpy array:

.. ipython:: python

    dip_traits_first_10 = np.array(pop.diploids.trait_array(np.array([i for i in range(10)])))

    dip_traits[:10] == dip_traits_first_10

You may also access with a slice:

.. ipython:: python

    dip_traits_first_10_via_slice = np.array(pop.diploids.trait_array(slice(0,10,2)))

    dip_traits_first_10[::2] == dip_traits_first_10_via_slice

The accesses are all range-checked and will raise exceptions when necessary:

.. ipython:: python

    #This tries to access past the end of the population:
    pop.diploids.trait_array(np.array([len(pop.diploids)+10]))

Slices are smart and can accomodate some types of input error:

.. ipython:: python

    #Oops--slice out of range!  But no exception will get raised:
    np.array(pop.diploids.trait_array(slice(len(pop.diploids)-10,len(pop.diploids)+10,3)))

Accessing diploid gamete keys
-------------------------------------------

You may also obtain structured arrays indicating the gametes in each individual:

.. ipython:: python

    gkeys = np.array(pop.diploids.key_array())
    gkeys.dtype

Just like the previous section, you may get data for subsets of the population using 
numpy arrays or slices.

Mutations
-------------------------------------------

You may get the mutations from the population into an array via:

.. ipython:: python

    muts = np.array(pop.mutations.array())
    muts.dtype

.. note::
    You may create a numpy array of the fixations list similarly.

The mutation keys in gametes
-------------------------------------------

You may access the mutation indexes in a gamete as follows:

.. ipython:: python

    #Get the neutral mutations 
    #from the 1st chromosome
    #of the 1st diploid:
    nkeys = np.array(pop.gametes[pop.diploids[0].first].mutations)
    nkeys.dtype

.. warning::
    Creating numpy arrays to the mutation keys gives you read-write
    access to the data!  This means that you can modify the state of the 
    population.  It is strongly recommended that you do not do this unless
    you know what you are doing :).  The following code block contains a way
    to prevent bad things from happening

.. ipython:: python

    nkeys.flags.writeable = False
    #Now, attempting to write raises an error
    nkeys[0] = 101

The mutation counts
-------------------------------------------

The mutation counts are stored using the same data type as the mutation keys in gametes:

.. ipython:: python

    mc = np.array(pop.mcounts)
    mc.dtype

.. note:: 
    The same caveat about read/write access applies here.

Example: sum of effect sizes in each gamete in a diploid
--------------------------------------------------------------------------------------

We'll integrate the above sections with an example calculating the sum of effect sizes
on each gamete in each diploid.  The result will be returned un a numpy structured array
that we pre-allocate to the correct size.

Here is our simple Python implementation.  We do not access the population data using
the buffer protocol.  The only use of numpy is to store the effect sizes for all mutations:

.. ipython:: python

    def get_esize_sum(gamete,esizes):
        s=0.0
        for i in gamete.smutations:
            s += esizes[i]
        return s

    def esize_sum_py(pop):
        rv=np.zeros((pop.N,2),dtype=[('first',np.float),('second',np.float)])
        i=0
        esizes=np.array([i.s for i in pop.mutations])
        for dip in pop.diploids:
            g1 = get_esize_sum(pop.gametes[dip.first],esizes)
            g2 = get_esize_sum(pop.gametes[dip.second],esizes)
            rv[i]=(g1,g2)
            i+=1
        return rv

The run times for the Python implementation:

.. ipython:: python

    %timeit esize_sum_py(pop)

Now, let's rewrite the above taking advantage of numpy arrays:

.. ipython:: python

    def get_esize_sum_np(keys,esizes):
        s=0.0
        for i in keys:
            s += esizes[i]
        return s

    def esize_sum_np(pop):
        rv=np.zeros((pop.N,2),dtype=[('first',np.float),('second',np.float)])
        d = np.array(pop.diploids.key_array())
        s = np.array(pop.mutations.array())['s']
        for i in range(d.shape[0]):
            g1 = get_esize_sum_np(np.array(pop.gametes[d['first'][i]].smutations,copy=False),s)
            g2 = get_esize_sum_np(np.array(pop.gametes[d['second'][i]].smutations,copy=False),s)
            rv[i]=(g1,g2)
        return rv

It is about twice as fast:

.. ipython:: python

    %timeit esize_sum_np(pop)

Check that both routines give the same answer:

.. ipython:: python

    check = esize_sum_py(pop)==esize_sum_np(pop)
    np.where(check == False)

The site-frequency spectrum
-------------------------------------------

This is the Python implementation from :ref:`processingpops`:

.. ipython:: python

    import collections

    def sfs_py(pop):
        sfsn = collections.Counter()   
        sfss = collections.Counter()   
        for i in range(len(pop.mcounts)):
            if pop.mcounts[i] > 0 :  
                if pop.mutations[i].neutral is False:
                    sfss[pop.mcounts[i]] += 1
                else:
                    sfsn[pop.mcounts[i]] += 1
        return (sfsn,sfss)

    %timeit sfs_py(pop)

Re-writing it to use numpy structured arrays is 10x faster.  The main 
difference is that we are using numpy to act as buffers to the underlying
C++ memory. Here it is:

.. ipython:: python

    def sfs_np(pop):
        m = np.array(pop.mutations.array())
        mc = np.array(pop.mcounts,copy=False)
        sfsn = collections.Counter()   
        sfss = collections.Counter()   
        extant = np.where(mc>0)
        n=m['neutral']
        for i in extant[0]: 
            if n[i]==1:
                sfsn[mc[i]] += 1
            else:
                sfss[mc[i]] += 1
        return (sfsn,sfss)

    %timeit sfs_np(pop)

They give the same results, too!

.. ipython:: python

    print(sfs_py(pop) == sfs_np(pop))

