.. _datamatrix:

Representing samples from a population in matrix form
====================================================================================

fwdpy11 allows viewing samples from populations (including the entire sample) in a matrix format.

The relevant functions are:

* :func:`fwdpy11.sampling.genotype_matrix`
* :func:`fwdpy11.sampling.haplotype_matrix`

Both return :class:`fwdpy11.sampling.DataMatrix` objects.

The API is extremely flexible, allowing for the inclusion or exclustion of selected markers.  Getting a 
:class:`fwdpy11.sampling.DataMatrix` requires a few steps:

1. Determine the individuals in your sample.  You may do this however you see fit.
2. Get a list of mutation keys associated with that sample via a call to :func:`fwdpy11.sampling.mutation_keys`.
3. Process/filter that list of keys as you need.  The list is unsorted.  You may choose to sort it according to mutation
   position, effect size, etc.  You may also apply a frequency filter.  The data returned from :func:`fwdpy11.sampling.mutation_keys` records the sample frequency and the `mcounts` objects present in populations can be used to filter on population frequencies.
4. Generate the fwdpy11.sampling.DataMatrix object by calling either :func:`fwdpy11.sampling.genotype_matrix` or
   `fwdpy11.sampling.haplotype_matrix`.
5. Get the matrices as numpy arrays via a call to :meth:`fwdpy11.sampling.DataMatrix.neutral` or :meth:`fwdpy11.sampling.DataMatrix.selected`.  These functions return 1d numpy arrays, which should be reshaped with the help of :attr:`fwdpy11.sampling.DataMatrix.nrow`, :meth:`fwdpy11.sampling.DataMatrix.ncol_neutral`, and :meth:`fwdpy11.sampling.DataMatrix.ncol_selected`.

The positions and population frequencies are also stored in the :class:`fwdpy11.sampling.DataMatrix` instance.  The
order of these mutations is the same order as the mutation keys used to generated the object.

The following example is a tour of the API:

.. testcode::

    import fwdpy11 as fp11
    import fwdpy11.wright_fisher as wf
    import fwdpy11.model_params
    import fwdpy11.sampling
    import numpy as np
    import pickle

    #First, we set up and run a 
    #simulation.
    N,theta,rho=1000,100,100

    p={'demography':np.array([N]*N,dtype=np.uint32),
       'nregions':[fp11.Region(0,1,1)],
       'recregions':[fp11.Region(0,1,1)],
       'sregions':[fp11.ExpS(0,1,1,0.25,0.25)],
       'rates':(theta/float(4*N),0.0,rho/float(4*N))
       }
    rng=fp11.GSLrng(42)
    params = fp11.model_params.SlocusParams(**p)
    pop=fp11.SlocusPop(N)
    #We simulate for N generations
    #because this code is run as part of the
    #testing suite, and so we want things
    #to be over quickly.
    pops = wf.evolve(rng, pop,params)

    #Now, we are going to represent the entire population
    #as a numpy matrix with dtype=np.int8.

    #Step 1.
    individuals=[i for i in range(pop.N)] #sample EVERYONE

    #Step 2.
    #By default, we get mutation keys back 
    #for neutral and selected mutations.
    #keys is a tuple.  keys[0] is neutral variants,
    #and keys[1] is selected variants
    keys = fp11.sampling.mutation_keys(pop,individuals)

    #Step3.
    #The keys come out totally unsorted.  Each element in
    #keys is itself a tuple.  The first element is the 
    #index of the mutation in pop.mutations and the 
    #second is the number of times it occurs in the sample
    #(which in this case is the entire population).
    #Let's sort the keys based on position and also remove singletons.
    neutral_sorted_keys=[i for i in sorted(keys[0],key=lambda x,m=pop.mutations: m[x[0]].pos) if i[1] > 1]
    selected_sorted_keys=[i for i in sorted(keys[1],key=lambda x,m=pop.mutations: m[x[0]].pos) if i[1] > 1]

    #Let's make sure we got that right:
    print(all(pop.mutations[neutral_sorted_keys[i][0]].pos <= 
        pop.mutations[neutral_sorted_keys[i+1][0]].pos for i in range(len(neutral_sorted_keys)-1)))
    print(all(pop.mutations[selected_sorted_keys[i][0]].pos <= 
        pop.mutations[selected_sorted_keys[i+1][0]].pos for i in range(len(selected_sorted_keys)-1)))

    #Step 4. -- get the DataMatrix encoded as a genotype matrix,
    #meaning 1 row per diploid and column values are 0,1,2
    #copies of derived allele
    dm = fwdpy11.sampling.genotype_matrix(pop,individuals,neutral_sorted_keys,selected_sorted_keys)

    print(type(dm))

    #Get the neutral genotypes out as a 2d 2d numpy array
    n = np.ndarray(dm.ndim_neutral(),buffer=dm.neutral,dtype=np.int8) 
    print(type(n))
    print(n.dtype)
    print(n.ndim)
    #This must be pop.N = 1,000:
    print(n.shape[0])

    #finally, the DataMatrix is picklable
    #As always with fwdpy11 types,
    #use -1 to select the latest
    #pickling protocol
    p = pickle.dumps(dm,-1)
    up = pickle.loads(p)

The output of the above code is:

.. testoutput::

    True
    True
    <class 'fwdpy11.sampling.DataMatrix'>
    <class 'numpy.ndarray'>
    int8
    2
    1000
