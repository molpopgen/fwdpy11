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

    rho=100.0

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
    assert len(x[0]) == 2*pop.N, "Node time error"

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


.. todo:: document limitations and future plans

Iterating over trees
------------------------------------

At the end of a simulation, a population's :class:`fwdpy11.ts.TableCollection` is 
population with a bunch of nodes, edges, etc..  But the "sequence" part of "tree
sequence" implies something about *iteration* that we haven't discussed yet.  fwdpy11
provides an efficient means of traversing the trees in a table collection in a 
left-to-right order along the genome.  The "visiting" of each tree is 
handled by :class:`fwdpy11.ts.TreeVisitor`, which gives you access to the 
:class:`fwdpy11.ts.MarginalTree` for each segment of the genome.

Traversing the trees is the core idea underying efficient algorithms for data analysis.
The multiply-linked list data structures stored in a :class:`fwdpy11.ts.MarginalTree` allow
for very rapid tree traversal.  Let's look at a concrete example.  We will calculate the 
average length of a marginal tree in our simulation.  To do this, we have to recognize the following: 

1. In these Wright-Fisher simulations, the tree times are in units of generations.
2. Each tree corresponds to a specific genomic segment, and these segment lengths differ
3. Thus, the mean total time on a tree is the weighted sum of the indiviudal marginal tree lengths.  
4. The weight on each tree is its genomic segment length divided by the genome length.
5. Here, the genome length is 1.0, which makes things easy (for once).

The numbers in comments at the ends of lines of code correspond to annotations following
immediately afterwards:

.. ipython:: python
    :okexcept:

    tv = fwdpy11.ts.TreeVisitor(pop.tables, [i for i in range(2*pop.N)])
    nodes = np.array(pop.tables.nodes, copy=False) # 1
    time = nodes['time'] # 1
    mean_total_time = 0.0
    while tv(False) is True: # 2
        m = tv.tree() # 3
        p = m.parents # 1
        segment_length = m.right - m.left
        tt_tree = 0.0
        for i in range(len(nodes)):
            if p[i] != fwdpy11.ts.NULL_NODE:
                branch_len = time[i] - time[p[i]] # 4
                mean_total_time += branch_len*segment_length
    print(mean_total_time/(4*pop.N))

1. We make several numpy arrays to view the data.  Internally, the data are stored in C++ containers.
   Thus, the numpy array is really a "view" of the data, and it requires no copies of the data.  However,
   It does take a small amount of time to make the view.  Thus, if we did not store the parents list in the 
   variable `p`, and instead referred to `m.parents` repeatedly, we would end up creating the view of the 
   parental data an additional `2*len(nodes)` times, and our calculation would slow down noticeably.
2. The `False` passed to the `__call__` function means "do not update the sample lists" for each tree.  The leaf
   count lists are always updated, however.  Saying `True` here updates the sample lists.  Sample list updating
   is relatively costly, which is why it is optional.
3. Internally, our TreeVisitor stores a C++ representation of a MarginalTree.  Here, through some C++ magic
   by the authors of pybind11, we are getting copy-free access to that stored data.
4. Time is measured from *past* to *present*. (This is a difference from msprime.)

The above loop is "Python fast", meaning that it is a pretty good mix of Python and numpy.  The main performance hits in
code like this are the looping and the round-trip from Python to numpy when accessing indexes in the numpy arrays.
These two performance bottlenecks have nothing to do with fwdpy11.  Rather, they are what we expect.  To do better, one
turns to the normal tricks, such as using Cython to move the operations entirely down to C.

It is now a good time to point out that total time counting is built-in because it is such a common operation:

.. ipython:: python

    # Need to construct a new visitor, as ours is all "iterated out"
    # from above
    tv = fwdpy11.ts.TreeVisitor(pop.tables, [i for i in range(2*pop.N)])
    mean_total_time = 0.0
    while tv(False) is True:
        m = tv.tree() # 3
        segment_length = m.right - m.left
        mean_total_time += segment_length*m.total_time(pop.tables.nodes)

    print(mean_total_time/(4*N))


The above loop is almost entirely composed of C++-side operations, and is thus extremely fast.

Constructing the TreeVisitor
+++++++++++++++++++++++++++++++++++++++++++++++++++++++

In the above example, the tree visitor was initialized using a :class:`fwdpy11.ts.TableCollection`
and a list of samples.  By setting the samples list equal to :math:`[0,2N)`, we are initializing with
respect to the last generation of the simulation.  Thus, the tree traversal will be updating the trees
for the entire population.  If you wish to iterate over the trees corresponding to a subset of the last generation,
simply create the approprate list, noting that the list may not contain redundant node ids.

A second method of initializing a TreeVisitor involves passing in two sample lists.  The intent here is that
the first list corresponds to the current generation ("alive nodes") and the latter to preserved nodes ("ancient
samples").  When passing in two lists, the tree iteration scheme tracks leaf counts separately for these two lists, 
via the fields :attr:`fwdpy11.ts.MarginalTree.leaf_counts` and :attr:`fwdpy11.ts.MarginalTree.preserved_leaf_counts`.  We'll see this in action
below.

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
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. ipython:: python

    import fwdpy11.tsrecorders
    pop = fwdpy11.SlocusPop(N,1.0)
    times = [5000]
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


Obtaining genotype data from tree sequences
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

You may obtain genotype data in the form of :class:`fwdpy11.sampling.DataMatrix`
objects. (See :ref:`datamatrix`.)  

The matrixes are constructed from a list of *node* ids (as opposed to individual indexes
as in :ref:`datamatrix`).  So, to get a matrix for our first 50 diploids:

.. ipython:: python

    m = fwdpy11.ts.make_data_matrix(pop, [i for i in range(100)], True, True)

The last two boolean arguments are whether or not to build data for neutral and selected
sites, respectively.  

The return value contains the following:

.. ipython:: python

    print(np.array(m.neutral))
    print(np.array(m.selected))

The neutral block is empty because we haven't added neutral mutations to our tree sequence yet (see above).

You may access the various fields using the usual operations:

.. ipython:: python

    print(m.selected.positions[:5])
    for i in range(5):
        assert m.selected.positions[i] == pop.mutations[m.selected_keys[i]].pos

The node id list passed to :func:`fwdpy11.ts.make_data_matrix` may contain nodes for alive samples 
or for ancient samples, allowing you to compare modern vs ancient nodes on the trees.

To get a list of node ids for ancient samples that are in the same order as the individuals to
which they belong, the following trick helps:

.. ipython:: python

    # Revisit our ancient node data
    print(ar[:5])
    # Stack and flatten gives us the nodes
    # for the individuals in the same order
    # as in ar:
    anodes = np.stack((ar['n1'],ar['n2']), axis=1).flatten()
    print(anodes[:10])


Tracking leaf counts separately for preserved and alive samples
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. todo:: show example

Initializing a simulation using msprime
------------------------------------------------------------------------

.. todo:: needs more testing

Outputting tables to msprime
------------------------------------------------------------------------

.. todo::

    Need to work out tables for individual metadata
