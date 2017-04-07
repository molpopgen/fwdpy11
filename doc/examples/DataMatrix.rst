Representing samples from a population in matrix form
====================================================================================

.. testcode::

    import fwdpy11 as fp11
    import fwdpy11.wright_fisher as wf
    import fwdpy11.sampling
    import numpy as np

    #First, we set up and run a 
    #simulation.
    N,theta,rho=1000,100,100

    nlist=np.array([N]*N,dtype=np.uint32)
    nregions=[fp11.Region(0,1,1)]
    recregions=nregions
    sregions=[fp11.ExpS(0,1,1,0.25,0.25)]
    rng=fp11.GSLrng(42)

    pop=fp11.Spop(N)
    #We simulate for N generations
    #because this code is run as part of the
    #testing suite, and so we want things
    #to be over quickly.
    pops = wf.evolve_regions(rng, pop,nlist, theta/float(4*N), 0, rho/float(4*N), nregions, [], recregions)

    #Now, we are going to represent the entire population
    #as a numpy matrix with dtype=np.int8.

    #Getting matrices out requires a few steps:
    #1. Determine the individuals in your sample.
    #2. Get a list of mutation keys associated with that sample.
    #3. Process/filter that list of keys as you need.
    #4. Generate the fwdpy11.sampling.DataMatrix object
    #5. Get the matrices ad numpy arrays

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

    #Get the neutral genotypes out as a numpy array
    n = dm.neutral()
    print(type(n))
    print(n.dtype)

    #n is a 1d array, and we want a 2d array
    #with rows as individuals and columns
    #as sites
    n = np.reshape(n,(dm.nrow,dm.ncol_neutral()))
    #This must be pop.N = 1,000:
    print(dm.nrow)

.. testoutput::

    True
    True
    <class 'fwdpy11.sampling.DataMatrix'>
    <class 'numpy.ndarray'>
    int8
    1000
