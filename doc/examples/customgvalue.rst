.. _customgvalues:

Custom stateless genetic value calculations
------------------------------------------------------------

Background:

* :ref:`definitions`

Further reading:

* For an example of a stateful model implemented in C++, see :ref:`stateful_fitness`.

fwdpy11 allows custom genetic value/fitness calculations to be implemented in either Python or C++.  The latter will
greatly outperform the former and we assume that Python implementations will be used for prototyping before writing
any C++ code.

The unit tests provide an example of implementing a fitness function where the fitnesses are multiplicative over
:math:`1`, :math:`1+t`, and :math:`1+s` for the three genotypes.  The trick we use is to realize that the `h` field of a
:class:`fwdpy11.Mutation` can be used as the `t`.  See the file `tests/test_custom_stateless_fitness.py` in the fwdpy11
repository for more details.

