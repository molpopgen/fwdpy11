.. _stateful_fitness:

Example of "stateful fitness": the snowdrift model
====================================================================

The low-level details in C++ are:

.. literalinclude:: ../../tests/ll_snowdrift.cc

The user-facing Python class is implemented using ``attrs``, which
is handy because we don't have to write C++ code to pickle/unpickle:

.. literalinclude:: ../../tests/snowdrift.py
