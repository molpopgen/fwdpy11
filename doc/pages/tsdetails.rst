.. _tsdetails:

Considerations to make when recording tree sequences
======================================================================

While tree sequence recording can be used to vastly speed up simulations, there are decisions that you can make 
that will speed things up even further.

The default behavior
-------------------------------------------------

What is happening behind the scenes of :func:`fwdpy11.wright_fisher_ts.evolve` when simplification occurs?
The default behavior is the following:

1. The tables are simplified and any node remapping of ancient samples occurs.
2. The new edge table is *indexed* according to the procedure described in Kelleher *et al.* (2016) PLoS Computational
   Biology, a.k.a the "msprime paper".
3. The number of occurrences of each mutation is counted by a tree traversal that follows the description of Algorithm L
   of the "msprime paper".
4. The mutation counts from step 3 are used to mark mutations for "recycling".  For example, a mutation with a count of
   zero has gone extinct.  Internally, fwdpp will re-use the memory of that extinct mutation.  Likewise, if fixations
   are being "pruned" from a simulation, fixed variants will be similarly recycled.

All of the above steps rely on efficient algorithms.  But not all of the above steps are *necessary*.

An alternative behavior
-------------------------------------------------

Typically, we will only simplify every 100 generations or so.  Thus, the mutation counts are only up to date with
respect to the last time simplification happened.  In this case, the mutations counts are simply a data structure being
used internally for orchestrating mutation recycling.  Related to this is that the table indexing is relatively
expensive, requiring two :math:`E\times log(E)` operations, where :math:`E` is the length of the edge table.

Steps two through four above may be skipped entirely by passing the value `True` to the `suppress_table_indexing`
argument of :func:`fwdpy11.wright_fisher_ts.evolve`.  This option triggers a code path that uses the return value from
simplification to mark mutations for recycling using a method that avoids the tree traversal entirely.

When should you use this option?  When some or all of the following apply:

1. You are simplifying very often but are not doing any time-series analysis during the simulation that relies on the
   mutation counts [1]_ [2]_.
2. You are running a simulation where you do not want to prune selected fixtations.  This option suppresses that
   possibility.
3. Your analysis of any genotypes will be restricted to the final generation of the simulation plus and ancient samples
   recorded during the simulation.

.. note::
    
    Suppressing table indexing is a huge performance boost when "densely" recording ancient samples. (For example,
    when taking a time series involving a random sample taken each generation.)  The reason is that ancient samples
    increase the size of the edge table because we must track the extra edges leading to those extinct lineages.  When
    tracking large numbers of ancient samples, the edge table grows considerably, meaning that traversal takes much
    longer, meaning that the default behavior described above is the less efficient option.


.. [1] If you insist, you may generate the mutation counts for the current generation by manually traversing all genomes
    for all diploids.

.. [2] Note that you may still analyze all individual metadata each generation.  Those data are always kept up to date.

Ending simulations early
-------------------------------------------------

.. versionadded:: 0.3.0

It is sometimes useful to simulate until a certain contition is met. fwdpy11 allows you to pass a callable to 
:func:`fwdpy11.wright_fisher_ts.evolve` that will return `True` if the desired condition is met and, if so, end the
simulation.  Let's consider a concrete example.  For example, when simulation a model of quantitative traits evolving to
an optimum shift, you may wish to end the simulation when the mean trait value gets above a certain value.  We can
represent this criterion via the following class:

.. code-block:: python

    import numpy as np

    class TraitValueHitsX(object):
        def __init__(self, threshold):
            self.threshold = threshold

        def __call__(population, simplified):
            md = np.array(population.diploid_metadata)
            if md['g'].mean() >= self.threshold:
                return True
            return False

The second parameter, `simplified` is `True` if the population's table collection has been simplified.
Stopping criteria depending on properties of mutations will only be efficient if this value is `True`.
Currently, such conditions are a bit tricky to handle, but one solution is to switch the simplification interval
to a smaller value when needed.

It is also possible to write these criteria as plain Python functions or as any C++ type convertible to
`std::function<bool(const fwdpy11::Population &, const bool)>`.
            
