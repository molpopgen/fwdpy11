Running simulations in parallel
==========================================

These scripts come from the fwdpy11 test suite.  They are integration tests to show that the various package components
are working together.  They also serve as nice examples.

Using Python's built-in methods for running multiple processes
-------------------------------------------------------------------------------

This first example shows the general recipe:

* Define a function to run your sims.  The function takes a `tuple` as arguments.
* Use Python 3's concurrent.futures_ to run multiple simulations in parallel

.. literalinclude:: ../../tests/test_simple_parallel.py


.. _concurrent.futures: https://docs.python.org/3/library/concurrent.futures.html

Recording the site-frequency spectrum every generation
-------------------------------------------------------------------------------

This is our function to evolve with a sampler:

.. literalinclude:: ../../examples/evolve_with_sampler.py

.. literalinclude:: ../../examples/async_sampler1.py
