.. _rng:

Random number generator
======================================================================

The random number generator type is :class:`fwdpy11.fwdpy11_types.GSLrng`.  It is initialized with a non-negative
integer as a seed:

.. testcode::

    import fwdpy11
    rng = fwdpy11.GSLrng(42)
