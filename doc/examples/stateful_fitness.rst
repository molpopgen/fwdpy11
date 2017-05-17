.. _stateful_fitness:

Custom stateful fitness models
==========================================

fwdpy11 allows custom fitness models to be implemented using C++.  Further, the fitness models may be stateful, meaning
that they keep track of their own data.  Shown below is an example of a "snowdrift" model.

To see an example of how the custom module is compiled and used, see the unit test file
`tests/test_stateful_fitness.py`.

.. literalinclude:: ../../tests/snowdrift.cpp

