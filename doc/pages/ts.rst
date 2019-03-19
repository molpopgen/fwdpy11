.. _ts:

Tree sequence recording
======================================================================

.. versionchanged:: 0.3.0

    The types :class:`fwdpy11.ts.TreeVisitor` and :class:`fwdpy11.ts.MarginalTree` are now deprecated, and their
    combined functionality are found in the new type, :class:`fwdpy11.ts.TreeIterator`, which provides much
    more efficient access to the marginal tree data.

Background reading:

1. :ref:`ts_data_types` gives links to the main data types involved.

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

    pop = fwdpy11.DiploidPopulation(N,1.0)

    rng = fwdpy11.GSLrng(42)

    rho=100.0

Set up a sim with no neutral regions.
Currently, attempting to simulate
neutral variants will throw an error because
I've not put some of the requisite tooling 
into the back-end, but that will come soon.

.. ipython:: python

    p = {'nregions':[],
    'sregions':[fwdpy11.GammaS(0,1,1,-10,1,scaling=2*N)],
    'recregions':[fwdpy11.Region(0,1,1)],
    'rates':(0.0,1e-3,rho/float(4*N)),
    'gvalue':fwdpy11.genetic_values.DiploidMult(2.0),
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


Saving tree sequences to files
-----------------------------------------------

Tree sequences are member data of populations.  Thus, they are serialized along with the population when it is pickled.
See :ref:`pickling_pops` for more details.

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
handled by :class:`fwdpy11.ts.TreeIterator`, which gives you access to the 
marginal tree data for each segment of the genome.

Traversing the trees is the core idea underying efficient algorithms for data analysis.
The multiply-linked list data structures stored in a :class:`fwdpy11.ts.TreeIterator` allow
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

    ti = fwdpy11.ts.TreeIterator(pop.tables, [i for i in range(2*pop.N)])
    nodes = np.array(pop.tables.nodes, copy=False) # 1
    time = nodes['time'] # 1
    mean_total_time = 0.0
    while ti() is True: # 2
        p = ti.parents # 1
        segment_length = ti.right - ti.left
        tt_tree = 0.0
        for i in range(len(nodes)):
            if p[i] != fwdpy11.ts.NULL_NODE:
                branch_len = time[i] - time[p[i]] # 4
                mean_total_time += branch_len*segment_length
    print(mean_total_time/(4*pop.N))

1. We make several numpy arrays to view the data.  Internally, the data are stored in C++ containers.
   Thus, the numpy array is really a "view" of the data, and it requires no copies of the data.  However,
   It does take a small amount of time to make the view.  Thus, if we did not store the parents list in the 
   variable `p`, and instead referred to `ti.parents` repeatedly, we would end up creating the view of the 
   parental data an additional `2*len(nodes)` times, and our calculation would slow down noticeably.
2. This is the call to "advance to the next tree".  When no more trees remain, `False` gets returned.
   By default, only "leaf counts" are updated as trees are advanced.  To update the sample lists for each tree,
   pass `True` as the last parameter to the `TreeIterator` constructor.
3. Time is measured from *past* to *present*. (This is a difference from tskit.)

The above loop is "Python fast", meaning that it is a pretty good mix of Python and numpy.  The main performance hits in
code like this are the looping and the round-trip from Python to numpy when accessing indexes in the numpy arrays.
These two performance bottlenecks have nothing to do with fwdpy11.  Rather, they are what we expect.  To do better, one
turns to the normal tricks, such as using Cython to move the operations entirely down to C.

It is now a good time to point out that total time counting is built-in because it is such a common operation:

.. ipython:: python

    # Need to construct a new visitor, as ours is all "iterated out"
    # from above
    ti = fwdpy11.ts.TreeIterator(pop.tables, [i for i in range(2*pop.N)])
    mean_total_time = 0.0
    while ti() is True:
        segment_length = ti.right - ti.left
        mean_total_time += segment_length*ti.total_time(pop.tables.nodes)

    print(mean_total_time/(4*N))


The above loop is almost entirely composed of C++-side operations, and is thus extremely fast.

Constructing the TreeIterator
+++++++++++++++++++++++++++++++++++++++++++++++++++++++

In the above example, the tree visitor was initialized using a :class:`fwdpy11.ts.TableCollection`
and a list of samples.  By setting the samples list equal to :math:`[0,2N)`, we are initializing with
respect to the last generation of the simulation.  Thus, the tree traversal will be updating the trees
for the entire population.  If you wish to iterate over the trees corresponding to a subset of the last generation,
simply create the approprate list, noting that the list may not contain redundant node ids.

A second method of initializing a TreeIterator involves passing in two sample lists.  The intent here is that
the first list corresponds to the current generation ("alive nodes") and the latter to preserved nodes ("ancient
samples").  When passing in two lists, the tree iteration scheme tracks leaf counts separately for these two lists, 
via the fields :attr:`fwdpy11.ts.TreeIterator.leaf_counts` and :attr:`fwdpy11.ts.TreeIterator.preserved_leaf_counts`.  We'll see this in action
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
    pop = fwdpy11.DiploidPopulation(N,1.0)
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

    print(pop.ancient_sample_metadata)
    for i in pop.ancient_sample_metadata[:5]:
        print(i.nodes)

We may view the same data using a numpy array.  These data are the same
format as for "alive" individuals.

.. ipython:: python

    ar = np.array(pop.ancient_sample_metadata, copy=False)
    print(ar.dtype)
    print(ar[:5])

Obtaining genotype data from tree sequences
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

You may obtain genotype data in the form of :class:`fwdpy11.sampling.DataMatrix`
objects. (See :ref:`datamatrix`.)  

The matrixes are constructed from a list of *node* ids (as opposed to individual indexes
as in :ref:`datamatrix`).  So, to get a matrix for our first 50 diploids:

.. ipython:: python

    m = fwdpy11.ts.data_matrix_from_tables(pop.tables, pop.mutations, [i for i in range(100)], True, True)

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

The node id list passed to :func:`fwdpy11.ts.data_matrix_from_tables` may contain nodes for alive samples 
or for ancient samples, allowing you to compare modern vs ancient nodes on the trees.

To get a list of node ids for ancient samples that are in the same order as the individuals to
which they belong, the following trick helps:

.. ipython:: python

    # Revisit our ancient node data
    print(ar[:5])
    anodes = ar['nodes'].flatten()
    print(anodes[:10])

You may also iterate over the genotypes on a per-marker basis using :class:`fwdpy11.ts.VariantIterator`:

.. ipython:: python

    vs = fwdpy11.ts.VariantIterator(pop)
    for v, i in zip(vs, pop.tables.mutations):
        r = v.record
        assert(v.genotypes.sum() == pop.mcounts[r.key])

The type of :attr:`fwdpy11.ts.VariantIterator.genotypes` is a numpy array with dtype `int8`.  The encoding is `0` for
the ancestral state and `1` for the derived.

Outputting tables to tskit
------------------------------------------------------------------------

You may convert the tables from a :class:`fwdpy11.ts.TableCollection` to a :class:`tskit.TreeSequence` as follows:

.. ipython:: python

    ts = pop.dump_tables_to_tskit()

The tables include information about mutations and individuals as metadata.  Each metadata record is a `UTF8`-encoded
string representations of `dict` objects.  Thus, to get the metadata back into something to work with, you may use
`eval`:

.. ipython:: python

    individial_zero_md = eval(ts.tables.individuals[0].metadata)
    print(individial_zero_md)

Note that it is only straightforward to get the first metadata record out!  The metadata columns are encoded as a big binary
blobs.  The tree sequence objects contain sets of vectors describing the *offset* of each metadata record, which is the
location of the start of each metadata record.  Thus, to access each record, you need to know its :math:`[start,end)`
position in the blob. Let's see how to use this and get the distribution of selection coefficients:

.. ipython:: python
    :okexcept:

    s = []
    for i in range(len(ts.tables.mutations)):
        j = ts.tables.mutations.metadata_offset[i]
        k = ts.tables.mutations.metadata_offset[i+1]
        d = eval(ts.tables.mutations.metadata[j:k])
        if d['neutral'] is False:
            s.append(d['s'])

    sa = np.array(s)
    print(sa.mean(), sa.min(), sa.max())

Let's sanity-check our result:

.. ipython:: python
    :okexcept:

    s2 = []
    for i in pop.tables.mutations:
        if pop.mutations[i.key].neutral is False:
            s2.append(pop.mutations[i.key].s)
    s2a = np.array(s2)
    print(s2a.mean(), s2a.min(), s2a.max())

Future releases of tskit will make metadata decoding a bit easier.

At this point, the major difference to be aware of is that the direction of time has been reversed.  With
that in mind, you may process the data in tskit, save it to a "trees file", etc.. See the tskit documentation_
for more details.

Tracking leaf counts separately for preserved and alive samples
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. todo:: show example

Initializing a simulation using msprime
------------------------------------------------------------------------

It is possible to start a simulation with a history generated by msprime_, via
:func:`fwdpy11.DiploidPopulation.create_from_tskit`

A large list of caveats apply:

1. Only nodes and edges are taken from an msprime TreeSequence.  If you wish to drop neutral mutation down, use
   :func:`fwdpy11.ts.infinite_sites`.
2. In general, it is difficult to know how long to simulate the resulting population, even when you start with
   a coalescent history.  Presumably, you will now simulate `pop` with selected variants and possibly with some
   "interesting" demography.  Be *very* careful about demography!  For example, let's assume that you simulate `pop` now
   with a fun demographic model for an additional `N` generations.  It is very unlikely that the population's ancestry
   will be independent of the original coalescent history.  What this means is that your intermediate-frequency variants 
   will be the product of a complex mixture of Kingman's coalescent process and the demographic model you applied.  In
   other words, your results will **not** reflect the equilibrium properties of the more complex model.
3. It is an error to simulate `pop` until all nodes originally present in `ts` are replaced and then stop the
   simulation, assuming you are now at equilibrium.  Such a procedure results in a biased sampling of the MRCA jump
   process, leaving you with tree sequences that are *shorter* than expected.
4. Probably more...

.. _documentation: https://msprime.readthedocs.io/en/stable/
.. _msprime: https://msprime.readthedocs.io/en/stable/
.. _tskit: https://tskit.readthedocs.io/en/stable/
