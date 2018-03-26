.. _sampling:

Sampling from simulated populations
======================================================================

To take a sample of :math:`n` diploids from a population, you have two ways to do things.  The first is to directly call
the `sample` member function of a population.  For example, :func:`fwdpy11.SlocusPop.sample`, as seen in the following
example:

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
    pdict = fwdpy11.ezparams.mslike(pop,simlen=pop.N, dfe=fwdpy11.ExpS(0,1,1,-0.1,1),pneutral = 0.95)
    params = fwdpy11.model_params.SlocusParams(**pdict)
    fwdpy11.wright_fisher.evolve(rng,pop,params)

Now, let's take a sample based on the first four diploids:

.. ipython:: python

    s = pop.sample(individuals = [0,1,2,3])

By default, s is a `tuple` with two elements.  The first represents neutral variants, and the second represents
selected. Let's look at the first few pieces of data for neutral variants:

.. ipython:: python

    print(s[0][:2])

Each variant is itself a tuple of (position, state).  The variant states are recorded as 0 = ancestral and 1 = derived.
Moving left to right along each string in each tuple, pairs of characters are the states corresponding to gametes 0 and
1 in individuals 0, 1, 2, and 3, respectively.  Across variant sites, the data are correctly phased.

These lists of tuples are compatible with pylibseq_ for processing.

To take a random sample *with* replacement from the population:

.. ipython:: python

    s2 = pop.sample(rng = rng, nsam = 4)

To take a random sample *without* replacement, numpy is your friend here:

.. ipython:: python

    np.random.seed(42)
    s3 = pop.sample(individuals = np.random.choice(pop.N, 4, replace=False))

If you do not want to separate neutral from selected variants:

.. ipython:: python

    s4 = pop.sample(individuals=[0,1,2,3], separate = False)

By default, mutations that are fixed in the sample are not included.  To change that behavior:

.. ipython:: python

    s5 = pop.sample(individuals=[0,1,2,3], remove_fixed = False)

You cannot mix and match the random sampling with a :class:`fwdpy11.GSLrng` and sampling a list of individuals:

.. ipython:: python

    s6 = pop.sample(rng = rng, individuals = [5,7,9,11])

Futher reading
-----------------------------------------------------

For finer-grained control over how samples are represented, see :ref:`datamatrix`.

.. _pylibseq: http://molpopgen.github.io/pylibseq/
