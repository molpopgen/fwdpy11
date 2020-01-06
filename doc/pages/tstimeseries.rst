.. _tstimeseries:

Efficient time series simulations with tree sequences.
====================================================================

Background reading:

* :ref:`timeseries`

Tree sequence recording provides a very efficient means of conducting a simulation and the tree sequence data structures are also very efficient for analyzing simulation results. In general, analyses using the tree
sequences will have logarithmic complexity (see Kelleher et al, 2016, PLoS Genetics).

When we are using ancient samples to "remember" or "preserve" large numbers of individuals during the simulation, the downstream processing is worse than logarithmic.  In the limiting case of preserving every individual ever simulated in a simulation wit ha constant population size, traversing through the trees has linear complexity.

In fwdpy11 0.5.2, we introduced the ability to "temporarily remember" individuals, meaning that we will "forget" them later in the simulation.
In essence, you may provide a callable function that processes the population immediately after simplification.
After execution of the function, all data concerning ancient samples will be cleared out from the table collection and from the population object.
This, these preserved nodes will be simplified out, or forgotten, the next time simplification occurs.

Let's look at a concrete example, where we track the frequency of all selected mutations over time.
First, we define a function to process our population and record our data.  We will track the unique
`(position, origin time, effect size)` tuples that uniquely identify a mutation, along with the frequency
associated with that tuple.  We will get these data by traversing all sample node time points, excluding
alive individuals:

.. ipython:: python

    from collections import namedtuple

    MutData = namedtuple('MutData',['pos','origin','s','daf'])

    freqs = []

    def track_freqs(pop):
       for t, n, m in pop.sample_timepoints(False):
          tables, idmap = fwdpy11.simplify_tables(pop.tables, n)
          trees = fwdpy11.TreeIterator(tables, idmap[n])
          for tree in trees:
             for mut in tree.mutations():
                k = mut.key
                p = pop.mutations[k].pos
                g = pop.mutations[k].g
                s = pop.mutations[k].s
                daf = tree.leaf_counts(mut.node)
                freqs.append(MutData(p,g,s,daf/len(n)))

The use of :func:`fwdpy11.DiploidPopulation.sample_timepoints` deserves some comment.  We pass `False`
to this function, telling it not to process the currently-alive individuals.  If we passed in `True` instead,
then we risk analyzing the nodes corresponding to those individuals twice, in the event that they get preserved
at some point in the future.

Now, we run a simple simulation with the above function passed into :func:`fwdpy11.evolvets`:

.. ipython:: python

    import fwdpy11
    import numpy as np

    tspop = fwdpy11.DiploidPopulation(100, 1.0)
    pdict = {'gvalue': fwdpy11.Multiplicative(2.),
              'nregions': [],
              'sregions': [fwdpy11.GammaS(0,1,1,2.0,1,scaling=200)],
              'recregions': [fwdpy11.PoissonInterval(0,1,0.1)],
              'rates': (0, 1e-2, None),
              'demography': np.array([100]*100, dtype=np.uint32)
              }
    params = fwdpy11.ModelParams(**pdict)

    preserver = fwdpy11.RandomAncientSamples(42, 100, np.arange(1, 100))

    rng = fwdpy11.GSLrng(125323)

    fwdpy11.evolvets(rng, tspop, params, 30,
                     recorder=preserver,
                     post_simplification_recorder=track_freqs,
                     suppress_table_indexing=True)

Finally, we plot our allele frequencies over time:

.. ipython:: python

    from matplotlib import rc
    rc('font',**{'size':18})
    import matplotlib.pyplot as plt
    import pandas as pd
    freqdf = pd.DataFrame(freqs, columns=MutData._fields)

    for n, g in freqdf.groupby(['pos','origin','s']):
        x = np.arange(len(g.daf))
        x += n[1]
        plt.plot(x, g.daf);

    plt.xlabel("Time (generation)");
    plt.ylabel("Mutation frequency");
    @savefig efficient_timeseries_example.png width=6in
    plt.tight_layout();
   
