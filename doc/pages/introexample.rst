Example
======================================================================

All of fwdpy11 is contained in a single import:

.. ipython:: python

    import fwdpy11
    import numpy as np

Let's set up some global variables for our simulation.  These are fairly standard:

.. ipython:: python

    N = 1000
    GENOME_LEN = 1.0
    THETA, RHO = 1000.0, 1000.0
    MU = 1e-3

We are simulating an additive effects model of a phenotype.  Fitness is determed by
the squared deviation from an optimum trait value using the classic model of Gaussian 
stabilizing selection.  In this simulation, the optimum will start out at zero and then change 
to 1 after `10N` generations of evolution.  We will use :class:`fwdpy11.GSSmo` to paramterize 
the optimum shift:

.. ipython:: python

    # The tuples are (generation, optimum, VS), where
    # VS is the inverse strength of stabilizing selection
    gssmo = fwdpy11.GSSmo([(0,0,1), (10*N,1,1)]) 


Our `gssmo` variable is an instance of a class that maps genetic values to fitness, and we use it to construct
our additive effect object, which is an instance of :class:`fwdpy11.Additive`:

.. ipython:: python

    additive_gss = fwdpy11.Additive(2.0, gssmo)

In the last cell, the `2` means that the model is additive over `0`, `sh`, and `2s` for `AA`, `Aa`, and `aa`
genotypes, respectively, where `a` is the mutant allele.  If you use a value of 1.0 instead, you will get the same
scaling as simulators like `slim`.  However, the 2.0 is useful here for the model of a quantitative trait.

.. note::

    The only difference between a simulation of a trait and a simulation of direct effects on 
    fitness is now the genetic value object is constructed. If no object derived from
    :class:`fwdpy11.GeneticValueToFitnessMap` is used, then the genetic value is treated as
    fitness itself.  In other words, if the previous cell had ommited the second argument, 
    then our simulation would be a standard simulation of additive effects on fitness:

    .. code-block:: python

        additive = fwdpy11.Additive(2.0)

Now, we use our parameters to construct an instance of :class:`fwdpy11.ModelParams`, which 
holds our parameters for us.  The `ModelParams` class takes `kwargs` as arguments. Our
preferred method for construction is to "explode" a `dict` containing our parameters:


.. ipython:: python

    pdict = {'gvalue': additive_gss,
            'nregions': [],
            'sregions': [fwdpy11.GaussianS(0, 1, 1, 0.15, 1)],
            'recregions': [fwdpy11.Region(0,1,1)],
            'rates': (0.0, MU, RHO/(4*N)),
            'demography': np.array([N]*(10*N + 100), dtype=np.uint32),
            'prune_selected': False
            }
    params = fwdpy11.ModelParams(**pdict)


Our population is an instance of :class:`fwdpy11.DiploidPopulation`:

.. ipython:: python

    pop = fwdpy11.DiploidPopulation(N, GENOME_LEN)

We also need a random number generator, which takes a 32-bit unsigned integer as a seed:

.. ipython:: python

    rng = fwdpy11.GSLrng(42)

fwdpy11 allows you to define arbitrary callables that process the population during simulation.
When recording tree sequences, a major use case for this processing is to define nodes to "preserve"
as "ancient samples".  What this means is that, at the end of the simulation, the nodes corresponding to 
these individuals will be retained in the tree sequences.  Their metadata will be preserved, too.

The callable must take two arguments. The first is the population, and the second is a Python object.  The 
type of the second argument's type is an internal detail.  It has a single user-facing interface, which is a function
called `assign`.  This function expects a numpy array (with a 32-bit signed integer dtype) containing the indexes of 
**individuals** to preserve.  Internally, these individual indexes will be converted to node indexes.

Below, we define a class that records **all** individuals in the population each generation after we have evolved to
equilibrium.  While we are at it, we will also record the generation and mean trait value, "because we can", and because 
it shows that we can basically do (almost) anything we want here in terms of time-series analysis.

.. ipython:: python

    class Recorder(object):
        def __init__(self, popsize):
            self.gbar = []
            self.individuals = np.arange(popsize, dtype=np.int32)
        def __call__(self, pop, ancient_sampler_recorder):
            if pop.generation >= 10*pop.N:
                md = np.array(pop.diploid_metadata, copy=False)
                self.gbar.append((pop.generation, md['g'].mean()))
                ancient_sampler_recorder.assign(self.individuals)

.. todo:: need to write a separate page on the details of recorders and tree sequences

.. ipython:: python

    recorder = Recorder(N)
    fwdpy11.evolvets(rng, pop, params, 100, recorder)

We can use the metadata to analyze our population. The metadata are represnted by 
the Python class :class:`fwdpy11.DiploidMetadata`.  The underlying C++ data 
structure is also registered as a numpy dtype, allowing more efficient analysis
through structured arrays.

.. ipython:: python

    # Let's get the mean trait value, the genetic variance and fitness
    # for the current generation
    alive_metadata = np.array(pop.diploid_metadata, copy=False)

    # The dtype names are the same as the DiploidMetadata 
    # class attributes.
    print(alive_metadata.dtype)

    # Note that alive_metadata does not own its data,
    # which means that the numpy array is just a
    # way to "view" to the C++ data without copying it.
    print(alive_metadata.flags)

Let's look at some properties of the final generation using both the Python class
and the structured array methods:

.. ipython:: python

    print(alive_metadata['g'].mean(), alive_metadata['g'].var(), alive_metadata['w'].mean())

    print(np.mean([i.g for i in pop.diploid_metadata]))
    print(np.var([i.g for i in pop.diploid_metadata]))
    print(np.mean([i.w for i in pop.diploid_metadata]))


Next, we will plot the mean trait value over time from the metadata.
The first thing we may want to take care of is that our metadata for 'alive'
and for 'ancient' samples are stored separately.  Let's fix that:

.. ipython:: python

    ancient_md = np.array(pop.ancient_sample_metadata, copy = False)
    all_md = np.concatenate((ancient_md, alive_metadata))

    # Combining the metadata resulted in a copy,
    # which you can see in the flags.  The new
    # object owns its data
    print(all_md.flags)

The access to fwdpy11 object data via numpy means that we can use the entire Python data stack.
Here, we will use `pandas` to get the mean trait value over time.  To do this, we first need 
the node times associated with our metadata nodes.  We will get these times by converting the population's
:class:`fwdpy11.NodeTable` into a structured array:

.. ipython:: python

    node_table = np.array(pop.tables.nodes, copy=False)
    print(node_table.dtype)
    mdtimes = node_table['time'][all_md['nodes'][:,0]]

Now, it is straightforward to create a `pandas.DataFrame` and aggregate with respect to time:
    
.. ipython:: python

    import pandas as pd
    df = pd.DataFrame(data={'time':mdtimes, 'g':all_md['g']})
    df = df.groupby(['time']).mean().reset_index()

The plotting is standard, too:

.. ipython:: python

    from matplotlib import rc
    rc('font',**{'size':18})
    rc('text', usetex=True)
    import matplotlib.pyplot as plt

    plt.plot(df.time, df.g);
    plt.ylabel("Mean trait value");
    plt.title("Adaptive walk to new optimum");
    plt.xlabel("Generation");
    @savefig mean_genetic_values_over_time.png width=6in
    plt.tight_layout();

Sanity check our calculations:

.. ipython:: python

    assert np.allclose(np.array([i[1] for i in recorder.gbar]), df.g) is True

.. ipython:: python

    ssh_over_time = []
    nmuts = fwdpy11.infinite_sites(rng, pop, THETA/(4*N))
    np.random.seed(54321)
    for t in np.unique(mdtimes):
        samples_at_t = np.where(mdtimes == t)[0]
        rsamples = np.random.choice(samples_at_t, 25, replace=False)
        rsamples_nodes = all_md['nodes'][rsamples,:].flatten()
        vi = fwdpy11.VariantIterator(pop.tables, pop.mutations, rsamples_nodes)
        ssh = 0.0
        for v in vi:
            g = v.genotypes
            r = v.record
            if pop.mutations[r.key].neutral is True:
                daf = float(g.sum())
                het = 2*daf*(len(g)-daf)/float(len(g)*(len(g)-1))
                ssh += het
        ssh_over_time.append(ssh)


    plt.plot(np.unique(mdtimes), ssh_over_time);
    plt.ylabel(r'$\pi$');
    @savefig pi_over_time.png width=6in
    plt.xlabel("Generation");
