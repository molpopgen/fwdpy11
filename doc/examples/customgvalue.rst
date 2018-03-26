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

fwdpy11's C++ headers provide macros to help you generate custom fitness functions.  There is support for fitness
calculations that depend on the genotype at a specific site as well as functions that act on a diploid as a whole.

The unit tests provide an example of implementing a fitness function where the fitnesses are multiplicative over
:math:`1`, :math:`1+t`, and :math:`1+s` for the three genotypes. This is an example of a fitness function acting on
genotypes at specific variable sites.  The trick we use is to realize that the `h` field of a
:class:`fwdpy11.Mutation` can be used as the `t`.  See the file `tests/test_custom_stateless_fitness.py` in the fwdpy11
repository for more details.  That file uses cppimport_ to compile a module based on the following C++ source code:

.. literalinclude:: ../../tests/custom_stateless_genotype.cpp

The unit tests also provide an implementation of a strict additive model without dominance.  The implementation acts on
the whole diploid. (It would be possible to have additivity with dominance, which would require the macros in the
previous example.) The unit test file is the same as above, but the C++ source is different:

.. literalinclude:: ../../tests/custom_additive.cpp

.. _cppimport: https://github.com/tbenthompson/cppimport
