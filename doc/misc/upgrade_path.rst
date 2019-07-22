.. _upgrade_path:

Upgrade path
====================================================================================

This document outlines how to upgrade existing scripts to new versions of fwdpy11.  This guide is likely
imperfect/incomplete.

0.5.0
-------------------------------------------------

The following functions and types previously required a :class:`fwdpy11.MutationVector` argument, but no longer do:

* :class:`fwdpy11.VariantIterator`
* :class:`fwdpy11.DataMatrixIterator`
* :func:`fwdpy11.data_matrix_from_tables`

The extra argument could be eliminated due to the new attributes added to :class:`fwdpy11.MutationRecord`.

0.2.0
--------------------------------------------------

This release also separates out the data representing a diploid into two classes, :class:`fwdpy11.DiploidGenotype` and
:class:`fwdpy11.DiploidMetadata`.  See :ref:`diploids` and :ref:`processingpopsNP` for type details and details on how
these new classes affect processing populations using NumPy, respectively.

This release contains major changes to how genetic values are calculated and to how simulations parameters are stored.
These changes are major *simplifications* to the package.  See :ref:`genetic_values_types` and :ref:`model_params` for
details.

The changes to how diploid data are stored completely changes how custom genetic values calculations are implemented.
See :ref:`customgvalues` and :ref:`stateful_fitness` for examples.

Another major change is that genetic value and noise functions are no longer allowed to be written in Python.  We may
bring that back in a later release.

class:`fwdpy11.sampling.DataMatrix` has been completely refactored.  See :ref:`datamatrix` for overview of current API.

The function :func:`fwdpy11.sampling.matrix_to_sample` now returns a tuple with two elements, which represent neutral
and selected gentoypes, respectively.  The previous  API made you choose neutral or selected for the return value, which
was a list.

Support for tree sequences will likely have a big impact on how you think about carrying out simulations.  See :ref:`ts`
and :ref:`ts_data_types` for details.

0.1.4
-----------------------------------

Changes to DataMatrix
+++++++++++++++++++++++++++++++++++++++

The member types :attr:`fwdpy11.sampling.DataMatrix.ndim_neutral` and  :attr:`fwdpy11.sampling.DataMatrix.ndim_selected` are now read-only attributes.  In previous versions, they were functions.  To upgrade, simply remove any trailing ``()``. In other words change this:

.. code-block:: python

   x.ndim_neutral()

To this:

.. code-block:: python

   x.ndim_neutral

The properties :attr:`fwdpy11.sampling.DataMatrix.neutral` and :attr:`fwdpy11.sampling.DataMatrix.selected` are now
writeable.  This allows you to recode the data as needed.  For example, if you wish to swap the 0/1 values for a column,
subtract 1 then multiply by -1.  The result will affect the data stored on the C++ side.


