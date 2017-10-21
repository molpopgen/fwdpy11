Upgrade path
====================================================================================

This document outlines how to upgrade existing scripts to new versions of fwdpy11.  This guide is likely
imperfect/incomplete.

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


