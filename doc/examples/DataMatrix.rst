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
5. Get the data as buffers via :attr:`fwdpy11.sampling.DataMatrix.neutral` or :attr:`fwdpy11.sampling.DataMatrix.selected`.  These buffers may be converted to numpy arrays with the help of :attr:`fwdpy11.sampling.DataMatrix.ndim_neutral`, :attr:`fwdpy11.sampling.DataMatrix.ndim_selected`.

The following is optional:

6. You may convert your resulting matrix into a sample format compatible with pylibseq_ via a call to
   :func:`fwdpy11.sampling.matrix_to_sample`.  Doing so only makes sense for a haplotype matrix, given how pylibseq_
   functions work.  For multilocus simulations, the output from :func:`fwdpy11.sampling.matrix_to_sample` may be
   separated out on a per-locus basis via a call to :func:`fwdpy11.sampling.separate_samples_by_loci`.

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

    # First, we set up and run a 
    # simulation.
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
    # We simulate for N generations
    # because this code is run as part of the
    # testing suite, and so we want things
    # to be over quickly.
    pops = wf.evolve(rng, pop,params)

    # Now, we are going to represent the entire population
    # as a numpy matrix with dtype=np.int8.

    # Step 1.
    individuals=[i for i in range(pop.N)] #sample EVERYONE

    # Step 2.
    # By default, we get mutation keys back 
    # for neutral and selected mutations.
    # keys is a tuple.  keys[0] is neutral variants,
    # and keys[1] is selected variants
    keys = fp11.sampling.mutation_keys(pop,individuals)

    # Step3.
    # The keys come out totally unsorted.  Each element in
    # keys is itself a tuple.  The first element is the 
    # index of the mutation in pop.mutations and the 
    # second is the number of times it occurs in the sample
    # (which in this case is the entire population).
    # Let's sort the keys based on position and also remove singletons.
    neutral_sorted_keys=[i for i in sorted(keys[0],key=lambda x,m=pop.mutations: m[x[0]].pos) if i[1] > 1]
    selected_sorted_keys=[i for i in sorted(keys[1],key=lambda x,m=pop.mutations: m[x[0]].pos) if i[1] > 1]

    # Let's make sure we got that right:
    print(all(pop.mutations[neutral_sorted_keys[i][0]].pos <= 
        pop.mutations[neutral_sorted_keys[i+1][0]].pos for i in range(len(neutral_sorted_keys)-1)))
    print(all(pop.mutations[selected_sorted_keys[i][0]].pos <= 
        pop.mutations[selected_sorted_keys[i+1][0]].pos for i in range(len(selected_sorted_keys)-1)))

    # Step 4. -- get the DataMatrix encoded as a genotype matrix,
    # meaning 1 row per diploid and column values are 0,1,2
    # copies of derived allele
    dm = fwdpy11.sampling.genotype_matrix(pop,individuals,neutral_sorted_keys,selected_sorted_keys)

    print(type(dm))

    # Get the neutral genotypes out as a 2d 2d numpy array
    n = np.ndarray(dm.ndim_neutral,buffer=dm.neutral,dtype=np.int8) 
    print(type(n))
    print(n.dtype)
    print(n.ndim)
    # This must be pop.N = 1,000:
    print(n.shape[0])

    # The DataMatrix is picklable
    # As always with fwdpy11 types,
    # use -1 to select the latest
    # pickling protocol
    p = pickle.dumps(dm,-1)
    up = pickle.loads(p)

    # We can also modify the data
    # in the array via Python's 
    # buffer protocol. Using
    # copy=False will give
    # us a buffer where modifications
    # will be passed on to the C++
    # side.
    
    # First, we'll copy
    # the existing view. 
    orig = n.copy()

    assert(n.shape == orig.shape)
    assert(np.array_equal(n, orig) == True)

    # We will swap all 0 and 2 encodings
    # in the data:
    for i in np.hsplit(n, n.shape[1]):
        i -= 2
        i *= -1

    # OK, let's prove that we've modified the C++
    # side.  We'll do that by making a new view,
    # and compare it to our copy:
    n2 = np.ndarray(dm.ndim_neutral,buffer=dm.neutral,dtype=np.int8) 
    assert(np.array_equal(n2, orig) == False)

    # Our new view is equivalent to our modified
    # view:
    assert(np.array_equal(n, n2) == True)

The output of the above code is:

.. testoutput::

    True
    True
    <class 'fwdpy11.sampling.DataMatrix'>
    <class 'numpy.ndarray'>
    int8
    2
    1000

Let's talk about what we did in this example.  We used the Python buffer protocol to view the genotypes at neutral
variants.  The ``buffer=`` argument to ``np.ndarray`` means that our NumPy array is a thin wrapper on top of memory
allocated in C++, giving us read-write access to the data.  The fact that we have write access allows our recoding of
the columns to be propagated to the C++ side.  Further, these thin wrappers give us very fast access to the underlying
data.

There are several use cases for recoding the data.  A DataMatrix is encoded by number of copies of the derived allele.
However, it may be useful to encode by number of copies of the minor allele, or the ``+`` allele when modeling a
quantitative trait.  For such cases, you can selectively recode the columns on a case-by-case basis.

It is possible to get a thin wrapper that is not writeable.  Doing so let's you have both fast access and safety. Let's revisit the above example:

.. ipython:: python
    :suppress:

    import fwdpy11 as fp11
    import fwdpy11.wright_fisher as wf
    import fwdpy11.model_params
    import fwdpy11.sampling
    import numpy as np
    import pickle

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
    wf.evolve(rng, pop, params)

    keys = fwdpy11.sampling.mutation_keys(pop, range(10))
    dm = fwdpy11.sampling.genotype_matrix(pop,range(10),keys[0],keys[1])

.. ipython:: python

    # Use a different syntax, to show that 
    # there are > 1 way to do things with
    # NumPy
    n = np.array(dm.neutral,copy=False,dtype=np.int8).reshape(dm.ndim_neutral)

Mark our new array as read-only:

.. ipython:: python

    n.flags.writeable = False

Now, we'll get an exception trying to modify the array:

.. ipython:: python

    for i in np.hsplit(n, n.shape[1]):
        i -= 2
        i *= -1


An array with ``flags.writeable=False`` can still be reshaped.  The flag simply prevents the raw-data from being
over-written.  The main use case for making an array read-only is to add a sense of safety to your code.  For example,
such arrays cannot be modified by functions.

.. _pylibseq: http://molpopgen.github.io/pylibseq/
