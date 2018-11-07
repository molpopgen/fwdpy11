.. _ts:

Tree sequence recording
======================================================================

Background reading:

1. :ref:`ts_data_types` gives links to the main data types involved.

.. todo:: Section on serialization

.. note::

    This page is currently limited in scope.  We'll add in more
    details in later releases.

Overview
------------------------------------

Example simulation
+++++++++++++++++++++++++++++++++++++++++

Let's set up a fast simulation, so that we have some data to work with in the later sections.

.. ipython:: python

    import fwdpy11
    import fwdpy11.model_params
    import fwdpy11.genetic_values
    import fwdpy11.ts
    import fwdpy11.wright_fisher_ts
    import numpy as np
    import pickle

We will simulate 1,000 diploids:

.. ipython:: python

    N = 1000

To initialize a population for tree sequence recording, pass a "genome length" 
as the second parameter to the Population constructor.This length is used to initialize
a :class:`fwdpy11.ts.TableCollection`

.. ipython:: python

    pop = fwdpy11.SlocusPop(N,1.0)

    rng = fwdpy11.GSLrng(42)

    rho=1000.0

Set up a sim w/no neutral parameters
Currently, attempting to simulate
neutral variant will throw an error because
I've not put some of the requisite tooling 
into the back-end, but that will come soon.

.. ipython:: python

    p = {'nregions':[],
    'sregions':[fwdpy11.GammaS(0,1,1,-10,1,scaling=2*N)],
    'recregions':[fwdpy11.Region(0,1,1)],
    'rates':(0.0,1e-3,rho/float(4*N)),
    'gvalue':fwdpy11.genetic_values.SlocusMult(2.0),
    'prune_selected': False,
    }
    params = fwdpy11.model_params.ModelParams(**p)
    params.demography=np.array([N]*10*N,dtype=np.uint32)

Run it, simplifying every 100 generations.  Note that we are calling 
:func:`fwdpy11.wright_fisher_ts.evolve`: to evolve the population.

.. ipython:: python

    fwdpy11.wright_fisher_ts.evolve(rng,pop,params,100)

Viewing a TableCollection
------------------------------------

Now, our object `pop` contains a tree sequence encoded in tables stored in a :class:`fwdpy11.ts.TableCollection`.  The
various tables can be viewed as numpy record arrays.  Viewing the tables this way gives access to the "bare" C++ types
via a "thin" wrapper in Python (the numpy array itself).

Let's take a look at the various tables.

.. ipython:: python

    # Don't forget the copy=False if you want MAX PERFORMANCE
    node_view = np.array(pop.tables.nodes, copy=False)
    print(node_view.dtype)
    print(node_view)

There must be 2N nodes with time equal to the current generation:

.. ipython:: python

    x = np.where(node_view['time'] == pop.generation)
    assert(len(x[0]) == 2*pop.N)

We may also look at the edges in the tree sequence:

.. ipython:: python

    edge_view = np.array(pop.tables.edges, copy=False)
    print(edge_view.dtype)
    print(edge_view)

We can get some useful info from the edges table.  For example, how many marginal trees are there?

.. ipython:: python

    print(len(np.unique(edge_view['left'])))

Finally, the mutation table:

.. ipython:: python

    mut_view = np.array(pop.tables.mutations, copy=False)
    print(mut_view.dtype)

The `key` field is the index of the mutation in the population's mutation vector:

.. ipython:: python

    for i in mut_view['key'][:5]:
        m = pop.mutations[i] 
        print(m.s,m.h,pop.mcounts[i])

Adding neutral mutations to a TableCollection
------------------------------------------------------------------------

So far, our population doesn't have any neutral variants.  Let's fix that:

.. ipython:: python

    theta = rho
    nmuts = fwdpy11.ts.infinite_sites(rng, pop, theta/(4*pop.N))
    print(nmuts)
    # have to recreate our view to the mutation table:
    mut_view = np.array(pop.tables.mutations, copy=False)
    for i in mut_view['key'][:5]:
        m = pop.mutations[i] 
        print(m.s,m.h,pop.mcounts[i])

We just added about 8,000 neutral mutations!

.. todo:: document limitations and future plans

Iterating over trees
------------------------------------

.. todo:: 

    Expose fwdpp's table_iterator to Python

Recording ancient samples during a simulation
------------------------------------------------------------------------

One of the selling points of tree sequences is a very efficient new method
of analyzing time series samples from simulations.  Without tree sequences,
we used "recorder" classes to analyze our populations during evolution. See
:ref:`recorders` for details.

Recording samples with tree sequences deciding which individuals to record
and then passing their indexes on to an instance of 
:class:`fwdpy11.tsrecorders.SampleRecorder`.  This class is a bridge between
you and the C++ back end.  The best way to show how to cross that bridge is 
to provide an example of a class that will take random samples from the
population at user-specified time points:

.. code-block:: python

    class RandomSamples(object):
        def __init__(self, nsam, timepoints):
            self.nsam = nsam
            self.timepoints = timepoints

        def __call__(self, pop, sr):
            if len(self.timepoints) > 0 and pop.generation > 0:
                if pop.generation == self.timepoints[0]:
                    # Make list of possible samples.  
                    # Note the dtype.
                    ind = np.arange(0, pop.N, dtype=np.uint32)
                    # requires that numpy be seeded
                    s = np.random.choice(ind, self.nsam, replace=False)
                    # Assign data to the SampleRecorder
                    sr.assign(s)
                self.timepoints.pop(0)

.. note::

    It is an error to attempt to preserve individuals from the final generation
    of a simulation as ancient samples.

The need to take random samples is so common that a class to do this is already provided.
See :class:`fwdpy11.tsrecorders.RandomAncientSamples` for details. We will use this built-in
type in the following section.

Viewing data for ancient samples
------------------------------------------------------------------------

.. ipython:: python

    import fwdpy11.tsrecorders
    pop = fwdpy11.SlocusPop(N,1.0)
    times = [i for i in range(2000, 10*pop.N, 2000)]
    # Parameters are: seed, sample size, time points:
    rec = fwdpy11.tsrecorders.RandomAncientSamples(14351, 50, times)
    fwdpy11.wright_fisher_ts.evolve(rng, pop, params, 100, rec)

At the end of the simulation, our population has a list of nodes corresponding to its
ancient samples:

.. ipython:: python

    print(pop.tables.preserved_nodes[:10])

Their node times must conform to what we expect:

.. ipython:: python

    print(np.unique([pop.tables.nodes[i].time for i in pop.tables.preserved_nodes]))

We also have *metadata* associated with our ancient samples.  For example, we 
have a mapping from their nodes to what individual they were:

.. ipython:: python

    print(pop.ancient_sample_records)
    for i in pop.ancient_sample_records[:5]:
        # time, node 1, node 2
        print(i.time, i.n1, i.n2)

We may view the same data using a numpy array:

.. ipython:: python

    ar = np.array(pop.ancient_sample_records, copy=False)
    print(ar.dtype)
    print(ar[:5])

The other form of metadata is the same as for alive individuals:

.. ipython:: python

    print(type(pop.ancient_sample_metadata[0]))
    md = np.array(pop.ancient_sample_metadata)
    print(md.dtype)
    print(md[:5])

Initializing a simulation using msprime
------------------------------------------------------------------------

Outputting tables to msprime
------------------------------------------------------------------------
