.. _conventions:

Documentation conventions
======================================================================

For convenience, importing fwdpy11 brings the following modules into the fwdpy11 namespace:

.. literalinclude:: ../../fwdpy11/__init__.py
    :lines: 20-

This allows you to write code like this:

.. doctest::

    import fwdpy11
    p = fwdpy11.SlocusPop(1000)

However, in order to refer to that class in the documentation, we must use its fully-specified name, which is
:class:`fwdpy11.fwdpy11_types.SlocusPop`.

