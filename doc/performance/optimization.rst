.. _optimization:

Optimizing performance
======================================================================

To keep your simulations as fast as possible, keep the following in mind:

* When processing data from a population, take advantage of being able to access the data as a numpy record array
  whenever possible.  See :ref:`processingpopsNP`.
* When writing your own fitness functions, do so in C++. See :ref:`customgvalues` and :ref:`stateful_fitness`.

Cython
--------------------------------

Cython_ is a static compiler that translates a Python-like grammar into C or C++.  It is possible that you may be able
to squeeze out some extra performance by writing some functions using Cython_.  Depending on what you are doing,
recorder objects (see :ref:`recorders`) may be made faster with Cython_, but the degree of improvement depends on the
details.  If your recorder is already calling a bunch of fwdpy11 functions, those are mostly in C++ already and will be
hard to outperform.  One of the unit test files (`tests/test_wright_fisher.py`) contains an example of a recorder
written in Cython_:

.. literalinclude:: ../../tests/MeanFitness.pyx

An important limitation is that you cannot define a Cython_ function to take the C++ representation of a fwdpy11 type as
an argument.  The reason is that Cython_ and pybind11_ (which fwdpy11 is based on) are different ways of reflecting C++
types in Python, and it does not seem straightforward to have the to systems interoperate.  See this issue_ for more
details. Note, however, that this is probably not a big deal.  The pybind11_-defined types in fwdpy11 are very efficient
Python classes with efficient access to member data.

One can of course use Cython_ to write interesting functions.  For example, here's one to calculate the mean effect size
of a mutation:

.. literalinclude:: ../../tests/mean_sel_coeff.pyx

You can see that function in action in `tests/test_numpy.py`.

.. _Cython: http://www.cython.org
.. _pybind11: https://pybind11.readthedocs.io/en/stable/
.. _issue: https://github.com/pybind/pybind11/issues/522
