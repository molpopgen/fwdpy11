.. _cpprecorder:

Example of time series analysis implemented in C++
=============================================================================

This page illustrates a few concepts:

1. fwdpy11 is extensible using both C++ and Python.
2. Often, the amount of C++ code to write is about as long as the equivalent Python
3. How to build plugins/extensions using cmake.

Using C++, we will write a callable that will record the mean trait value 
of a population over time.  In essence, we will build a small Python 
extension module using pybind11 that will be called `gvalue_recorder`. The
module will define a single Python function, `gvalue_recorder.record_gvalue`.
That function will take a single Python list as an argument, and returns a Python
object known as a capsule_, which is a Python object holding something defined in 
another language (C++ in thise case).  That capsule holds our callable, which follows the requirements described in
:ref:`timeseries`.

Let's just look at the code:

.. literalinclude:: ../../examples/plugin/gvalue_recorder.cc

The code is very simple:

* The Python function is defined using a C++ lambda
* The Python function returns a `pybind11::cpp_function` (which is a capsule_) that contains another C++ lambda that
  "captures" the argument passed to the Python function.

We have a function returning a function, in other words.

The following Python script tests that our plugin behaves as expected:

.. literalinclude:: ../../examples/plugin/test_plugin.py

For the record, the Python equivalent to what we've generated in C++ is:

.. code-block:: python

    gvalues = []
    def r(pop, sampler):
        md = np.array(pop.diploid_metadata, copy=False)
        gvalues.append(md['g'].mean())

I expect the performance of both methods to be nearly equivalent, largely because traversing the metadata once per
generation is extremely fast compared to everything else going on in a simulation.

The plugin is built using cmake.  Note the `execute_process` steps, which find the fwdpy11 headers for you. 

.. literalinclude:: ../../examples/plugin/CMakeLists.txt

.. _capsule: https://docs.python.org/3/c-api/capsule.html
