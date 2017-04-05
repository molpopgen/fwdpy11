fwdpy11
*************************

This is the README for fwdpy11_, which is a Python package for forward-time population genetic simulation.  It uses
fwdpp_ as its C++ back-end.

Build status
-----------------------

.. image:: https://travis-ci.org/molpopgen/fwdpy11.svg?branch=master
    :target: https://travis-ci.org/molpopgen/fwdpy11
	:alt: Master branch

.. image:: https://travis-ci.org/molpopgen/fwdpy11.svg?branch=dev
    :target: https://travis-ci.org/molpopgen/fwdpy11
	:alt: Development branch

Manual
-----------------------

.. image:: https://readthedocs.org/projects/fwdpy11/badge/?version=latest
	:target: http://fwdpy11.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status (master branch)

Features
-----------------------

* Picklable population objects
* Parallel computation via multiprocessing_ or concurrent.futures_.
* Custom temporal samplers to analyze populations *during* a simulation may be written in pure Python.
  
Dependencies
-----------------------

The following must be present on your system:

* GSL_. This is a C library.  It is available via `conda`.
* pybind11_. This should be installed via `pip` or `conda` as appropriate for your system.

fwdpy11_ also uses fwdpp_, which is included as a submodule.

License
-----------------------

GPLv3 or later (See COPYING)

Suppored Python version
-----------------------

Python 3.5 or newer is preferred.  Most features will work with 2.7 or later, but you will not be able to retrieve
pickled populations from files using Python2.7.  

.. code-block:: bash

    git submodule init
    git submodule update
    python setup.py build_ext -i
    python test.py

.. _fwdpy11: http://molpopgen.github.io/fwdpy11
.. _fwdpp: http://molpopgen.github.io/fwdpp
.. _GSL: http://gnu.org/software/gsl
.. _pybind11: https://github.com/pybind/pybind11
.. _multiprocessing: https://docs.python.org/3/library/multiprocessing.html
.. _concurrent.futures: https://docs.python.org/3/library/concurrent.futures.html
