.. _sampling:

Sampling from simulated populations
======================================================================

.. versionchanged:: 0.2.0

    The return value of sampling routines changed to a :class:`fwdpy11.sampling.DataMatrix`

To take a sample of :math:`n` diploids from a population, you have two ways to do things.  The first is to directly call
the `sample` member function of a population.  For example, :func:`fwdpy11.SlocusPop.sample`, as seen in the following
example:

.. ipython:: python

    import fwdpy11
    import fwdpy11.genetic_values
    import fwdpy11.model_params
    import fwdpy11.ezparams
    import fwdpy11.wright_fisher
    import numpy as np
    rng = fwdpy11.GSLrng(42)
    pop = fwdpy11.SlocusPop(1000)
    pdict = fwdpy11.ezparams.mslike(pop,simlen=pop.N, dfe=fwdpy11.ExpS(0,1,1,-0.1,1),pneutral = 0.95)
    params = fwdpy11.model_params.ModelParams(**pdict)
    fwdpy11.wright_fisher.evolve(rng,pop,params)

Now, let's take a sample based on the first four diploids:

.. ipython:: python

    s = pop.sample(individuals = [0,1,2,3])
    print(type(s))

The return value is a :class:`fwdpy11.sampling.DataMatrix`. See :ref:`datamatrix` for more details.

By default, a *haplotype* matrix is returned:

.. ipython:: python

    hm = np.array(s.neutral)
    print (np.unique(hm))
    print (hm.shape)

To get a genotype matrix:

.. ipython:: python

    s = pop.sample(individuals = [0,1,2,3], haplotype = False)
    gm = np.array(s.neutral)
    print (np.unique(gm))
    print (gm.shape)

The :class:`fwdpy11.sampling.DataMatrix` format is compatible with scikit-allel_ and the new data structures that will
be found in libsequence_ 2.0 (and therefore a future version of pylibseq_).  To work with the current pylibseq_, you
need to convert the matrix format to a list of tuples:

.. ipython:: python

    # Need a haplotype matrix
    s = pop.sample(individuals = [0,1,2,3])
    neutral, selected = fwdpy11.sampling.matrix_to_sample(s)
    print(neutral[0])

The format of each tuple element is `(position, genotypes)`, where the genotypes are 0 = ancestral, 1 = derived, and the
first pair of values corresponds to the first individual, etc.

To combine the neutral and selected data into a single block usable in pylibseq_:

.. ipython:: python

    combined = neutral + selected
    # sort on positions
    combined.sort(key = lambda x: x[0])

.. _pylibseq: http://molpopgen.github.io/pylibseq/
.. _libsequence: http://molpopgen.github.io/libsequence
.. _scikit-allel: https://scikit-allel.readthedocs.io/en/latest/
