.. _genotypes_trees:


Accessing genotypes and individual trees
======================================================================

Iterating over individual variants
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Let's get the derived mutation frequency spectrum at neutral and selected sites
for the first fifty nodes in the population generated in :ref:`introexample`.
To do so, we use :class:`fwdpy11.VariantIterator` to traverse each mutation:

.. ipython:: python

   neutral_sfs = np.zeros(50)
   selected_sfs = np.zeros(50)
   vi = fwdpy11.VariantIterator(pop.tables, [i for i in range(50)])
   for v in vi:
       n = pop.mutations[v.records[0].key].neutral
       g = v.genotypes
       daf = g.sum()
       if n is True:
           neutral_sfs[daf - 1] += 1
       else:
           selected_sfs[daf - 1] += 1

Notes:

* `v.record` is the :class:`fwdpy11.MutationRecord` corresponding to the variant
* `g` is a numpy array of 8-bit integers. 0 = ancestral state and 1 = derived state. The element order 
  is the same as the order of sample indexes used to create `vi`.
* The iteration includes sites that are fixed in the sample.
* You may iterate over the variants in the genomic window :math:`[begin, end)` by passing the appropriate
  values to a :class:`fwdpy11.VariantIterator` constructor.

.. ipython:: python

   f, (ax1) = plt.subplots(1, 1)
   ax1.plot([i+1 for i in range(50)], selected_sfs/selected_sfs.sum(),
            label='selected');
   ax1.plot([i+1 for i in range(50)], neutral_sfs/neutral_sfs.sum(),
            label='neutral');
   ax1.set_xlabel("Derived allele frequency");
   ax1.set_ylabel("Proportion of variants");
   ax1.set_xlim(1,50);
   ax1.set_xticks([1,10,20,30,40,50]);
   plt.legend(loc='best');
   @savefig sfs_example.png width=6in
   plt.tight_layout()

Obtaining genotype data as a numpy matrix
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

We may obtain the genotypes for the samples all at once using the following function:

.. autofunction:: fwdpy11.data_matrix_from_tables

.. ipython:: python

   dm = fwdpy11.data_matrix_from_tables(pop.tables, [i for i in range(50)],
                                        record_neutral=False,
                                        record_selected=True,
                                        include_fixations=True)
   selected_sfs_from_matrix = np.zeros(50)
   # Get the row sums, which are the DAF
   # for each site
   rc = np.sum(dm.selected, axis=1)
   n, c = np.unique(rc, return_counts=True)
   for i, j in zip(n, c):
       selected_sfs_from_matrix[i-1] = j
   assert np.array_equal(selected_sfs, selected_sfs_from_matrix)

Notes:

* `dm` is a :class:`fwdpy11.DataMatrix`.
* By default, fixations in the sample are *not* included (which saves memory).  
  We included them here to test results against what we did above.
* You may iterate over the variants in the genomic window :math:`[begin, end)` by passing the appropriate
  values to :func:`fwdpy11.data_matrix_from_tables`.
* We skipped neutral mutations in this example to save memory on the servers that build this manual.

Accessing genotypes from multiple intervals via repeated calls to :func:`fwdpy11.data_matrix_from_tables` is an
:math:`O(n^2)` algorithm.  The reason is that traversal through the tree sequences will start over with each new call.
If you need to access genotypes from multiple genomic windows, you may use :class:`fwdpy11.DataMatrixIterator` instead, 
which will be much more efficient.  The following example extracts two windows corresponding to the same samples as 
the previous example, and tests that the data are identical:

.. ipython:: python

    windows = [(0, 0.25), (0.75, 1.0)]
    selected_genotypes = np.array(dm.selected)
    p = np.array(dm.selected.positions)
    dmi = fwdpy11.DataMatrixIterator(pop.tables,
                                    [i for i in range(50)], windows,
                                    neutral=False, selected=True, fixations=True)
    for i,j in zip(dmi, windows):
        idx = np.where((p >= j[0]) & (p < j[1]))[0]
        slice = selected_genotypes[idx,:] 
        assert np.array_equal(i.selected, slice)

Tree traversal
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

:class:`fwdpy11.TreeIterator` allows left-to-right traversal of the marginal trees:

.. ipython:: python

   ti = fwdpy11.TreeIterator(pop.tables, [i for i in range(50)])

For each tree, you may access the parent/child/etc. elements described in :ref:`representing_trees`.  In general,
tree access is an advanced topic that allows efficient algorithms to be developed.  For now, see :ref:`bgs` for an
example.

For each marginal tree, you have access to the following:

* The list of nodes, output in a preorder traversal via :func:`fwdpy11.TreeIterator.nodes`.
* The complete list of samples via :func:`fwdpy11.TreeIterator.samples`
* The list of samples below any node via :func:`fwdpy11.TreeIterator.samples_below`.
* Iterators to the sites and mutations on the current tree via :func:`fwdpy11.TreeIterator.sites` and :func:`fwdpy11.TreeIterator.mutations`, respectively.
