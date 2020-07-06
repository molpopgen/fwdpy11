.. _advancedtopics:

Advanced topics
=================================================

Executing multiple replicates in parallel using `concurrent.futures`
---------------------------------------------------------------------

You have many options when it comes to running many replicates of a model.
For example, you could write a single Python script plus a shell script
that your local cluster understands, sending your simulation to hundreds
of compute nodes.

This section introduces :mod:`concurrent.futures` as a method to execute
many replicates in separate ``Python`` *processes* using a single script.

Fully worked-out examples of varying complexity can be found here:

* :ref:`bgs`
* :ref:`IMexample`
* :ref:`migtest`
* :ref:`precapitation`
* :ref:`recapitation`

A common set of idioms used by these examples are:

* The function to run the simulation typically does not return
  the instance of :class:`fwdpy11.DiploidPopulation`.  Such
  a return requires pickling, which is very expensive.  If you
  need to store these population, write them to files rather
  than returning them. See :ref:`here <savingsimstodisk>`.
* If you do return something to the main "collector"
  process, keep it simple!
* Random number seeds are set in ``__main__`` using 
  :func:`numpy.random.randint` after calling
  :func:`numpy.random.seed` using a user-provided seed.

.. _processing_metadata:

Decoding metadata from :class:`tskit.TreeSequence`
---------------------------------------------------------------------

`tskit` and `fwdpy11` treat metadata quite differently.  The former is much more general,
while the latter gives you direct access to the data objects on the C++ side.
The `tskit` approach is based on binary strings.  What `fwdpy11` does is encode strings that
can be converted back to Python dictionaries.  For example, here is how one may process the
individual metadata after dumping the tables to `tskit`:

.. code-block:: python

   import tskit

   individuals = ts.tables.individuals
   md = tskit.unpack_bytes(individuals.metadata, individuals.metadata_offset)
   first_individual = eval(md[0])

In order to distinguish "alive" from "dead" individuals (*e.g.*, those 
preserved as ancient samples), we need to make use of flags found in
:mod:`fwdpy11.tskit_tools`. For example, to identify all currently alive
individuals in the individual table:

.. code-block:: python

    import fwdpy11.tskit_tools

    alive_individuals = (
        ts.tables.individuals.flags & fwdpy11.tskit_tools.INDIVIDUAL_IS_ALIVE
    )

To find the nodes corresponding to those individuals:

.. code-block:: python

    alive_individual_nodes = (ts.tables.nodes.individual >= 0) & (
        (
            ts.tables.individuals.flags[ts.tables.nodes.individual]
            & fwdpy11.tskit_tools.INDIVIDUAL_IS_ALIVE
        )
        > 0
    )

    # The time of these nodes should all be zero because we simulated
    # non-overlapping generations
    alive_node_times = np.unique(
        ts.tables.nodes.time[np.where(alive_individual_nodes)], return_counts=True
    )

.. todo::

    Provide nice functions to return nodes from various times,
    etc..

The mutation metadata follow the same general recipe:

.. code-block:: python

   mutations = ts.tables.mutations
   sites = ts.tables.sites
   md = tskit.unpack_bytes(mutations.metadata, mutations.metadata_offset)

Here, ``md`` is a :class:`dict` whose key names are the same as the attributes of 
:class`fwdpy11.Mutation`, with one exception. The `tskit` representation
of the mutation record's the allele's *age* in generations while
:attr:`fwdpy11.Mutation.g` is the generation when the mutation arose.
The reason for this discrepancy is because ``fwdpy11`` thinks forward in time
while ``tskit`` thinks backwards. The conversion to and from is trivial:

.. code-block:: python

   print(f"Origin to age = {pop.generation - m.g + 1}")

.. _howlongtorun:

How long to run the simulation
---------------------------------------------------------------------

Deciding how long to run a simulation depends on whether or not
a "burn in" phase is required to get the model to statistical
equilibrium.  For example, if you want to get the steady-state
properties of a model with many selected mutations and linkage,
you will need to burn in.  If you want to know how the short term
dynamics of the neutral site-frequency spectrum are affected by suddenly
introducing strongly-beneficial mutations, then you may not need to 
burn in.

So, then, how long to burn in?  Clearly, the answer should be "the 
least that you have to"!  As always, there are a few different
things to consider:

* You want the *distribution* of final outcomes to be "correct".
  Theory often gives us expressions for expectations and 
  sometimes for variances of things we care about, like
  the number of mutations in a sample. However, comparing
  distributions will often require a statistical analysis
  comparing the output of independent simulations.
  For example, for a neutral model, the distribution of the
  number of mutations under an infinitely-many sites model
  for a small sample should match ``msprime`` very well.
  For a steady-state model of recurrent hitch-hiking
  [KaplanHudsonLangley1989]_, results from small samples
  should match coalescent simulations ([KernSchrider2016]_
  for example). Clearly, all of the different methods
  should match analytical results where possible.
* You probably want all of your trees to be fully coalesced
  to a single common ancestor.  The time back to a final
  ancestor has a very long tail.  In models with recombination,
  it is not uncommon to have a tree or three with multiple roots
  at the end of a simulation.

For the first point, a burn-in of ``10N`` generations or so
is usually sufficient. Before we simulated with tree sequence 
recording, we simulated entire "genomes" (neutral mutations and all),
which was rather slow.  (Here "we" is the field in general.)  We then
(hopefully!) compared our results to something like ``msprime``.  The 
agreement was really good, so one answer is "about ``10N``.  I'm being
vague about ``N`` here--typically, it would be equal to the starting effective
population size of your population (*e.g.* in the absence of selection).

With tree sequence recording, we run into the "uncoalesced marginal
trees" issue mentioned above.  We could just simulate much longer
to get rid of this problem, but it is simply easier to start
with a tree sequence from ``msprime`` (see :ref:`here <starting_from_msprime>`).


.. _eyre_walker:

The "Eyre-Walker" model of complex traits
---------------------------------------------------------------------

[EyreWalker2010]_ describes a model relating mutations with
fitness effect :math:`S = |2Ns|` to :math:`z`, their effect on
a phenotype/trait according to :math:`z = \delta S^\tau(1+\epsilon)`.
Here, :math:`\epsilon` is a draw from a Gaussian distribution with mean zero,
:math:`\delta` is :math:`1` or :math:`-1` with equal probability, and 
:math:`\tau` "tunes" the correlation between fitness effect and effect on the trait.

Implementing this model is quite straightforward, as the trait values do not affect
the dynamics of the model.  For this flavor of a "pure pleiotropy" model,
the technical details reduce to a standard population genetic model
where we tack on the trait values at the end.

Let's set up a model where :math:`S` is exponentially distributed with mean ``-20``.
We'll run a small population for a few generations.  This model will be nowhere
near equilibrium, but we're just using it as an example:

.. ipython:: python

    import fwdpy11
    import numpy as np

    N = 1000
    sregions = [fwdpy11.ExpS(beg=0.0, end=1.0, weight=1.0, mean=-20, scaling=2 * N)]
    recregions = [fwdpy11.PoissonInterval(beg=0.0, end=0.1, mean=1e-3)]
    gvalue = fwdpy11.Multiplicative(scaling=2.0)
    pdict = {
        "nregions": [],
        "sregions": sregions,
        "recregions": recregions,
        "rates": (0.0, 1e-3, None),
        "gvalue": gvalue,
        "prune_selected": False,
        "simlen": 150,
    }
    params = fwdpy11.ModelParams(**pdict)
    pop = fwdpy11.DiploidPopulation(N, 1.0)
    rng = fwdpy11.GSLrng(54321)
    fwdpy11.evolvets(rng, pop, params, 100)
    print(len(pop.tables.mutations))


Define a function relating :math:`S` to :math:`z`:

.. ipython:: python

    def getz(S, tau, sigma):
        if np.random.uniform() < 0.5:
            delta = 1
        else:
            delta = -1
        epsilon = np.random.normal(loc=0, scale=sigma, size=1)[0]
        return delta * np.power(np.abs(S), tau) * (1.0 + epsilon)

Apply the function and look at the results:

.. ipython:: python

    np.random.seed(101010)
    zvals = {}
    for i, m in enumerate(pop.tables.mutations):
        zvals[m.key] = getz(2 * N * pop.mutations[m.key].s, 0.5, 1.0)

    for k, v in zvals.items():
        print(pop.mutations[k].s, v)

Traverse the tree sequence to get individual phenotypes under a 
strictly additive model:

.. ipython:: python

    phenotypes = np.zeros(pop.N)
    node_to_individual = {}
    for i, j in enumerate(pop.diploid_metadata):
        assert j.nodes[0] not in node_to_individual
        assert j.nodes[1] not in node_to_individual
        node_to_individual[j.nodes[0]] = i
        node_to_individual[j.nodes[1]] = i
    ti = fwdpy11.TreeIterator(pop.tables, pop.alive_nodes, update_samples=True)
    for t in ti:
        for m in t.mutations():
            for n in t.samples_below(m.node):
                phenotypes[node_to_individual[n]] += zvals[m.key]

The trait value distribution is:

.. ipython:: python

    np.unique(phenotypes, return_counts=True)

The mean trait value and the genetic variance are:

.. ipython:: python

    phenotypes.mean()
    phenotypes.var()

For our final trick, let's store them in the individual metadata:

.. ipython:: python

    md = np.array(pop.diploid_metadata, copy=False)
    md['g'][:] = phenotypes

    for i, j in zip(md['g'][:5], phenotypes[:5]):
        print(i, j)

The trick is that the ``numpy`` array is a *non-owning* array, meaning
that it is a simple "view" of the underlying ``C++`` data.  Further,
it happens to be read/write, and thus we can modify it.  (Try not
to abuse this.  More often than not, you'll just break stuff.)

I feel that some comments about this model are warranted:

* The model is a simplification of earlier work by Keightley and Hill
  ([HillKeightley1988]_, [KeightleyHill1990]_), which should be cited
  alongside [EyreWalker2010]_.
* [HillKeightley1988]_ and [KeightleyHill1990]_ show that this model predicts heritability 
  (the genetic variance) is linear-ish with ``N``, with the details depending
  somewhat on the model parameters.  The biological reasonableness
  of that prediction is dubious.  [JohnsonBarton2005]_ discuss this point
  in some detail, and give other relevant references.
* Depending somewhat on the parameters of the ``getz`` function, one will eventually
  generate an intermediate-frequency variant with a massive :math:`z`,
  meaning that it would explain a considerable amount of the genic
  variance for the trait (high :math:`2pqz^2`).  Such outcomes are contrary
  to the results of human GWAS.
* Many applications of this model in the literature do artitrary things to
  the simulations, like only treating mutations in certain frequency ranges
  as affecting trait values.  Thus, the predictions made by such studies are
  not natural outcomes of evolutionary models.  Consider only applying
  the ``getz`` functtion to mutations with frequency :math:`< x` in order
  to say something about "rare alleles".  This treatment of the data actually
  changes the model: mutations that are common now (at the end of the simulation)
  were rare at some point in the past, and had frequencies :math:`< x`.  Therefore,
  the model is one where mutations suddenly stop affecting the trait once
  they hit some critical frequency. Presumably, the same variants would affect
  the trait again should they drift to a frequency below :math:`x` if the
  simulation is run a bit longer.
