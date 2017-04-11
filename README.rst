fwdpy11
*************************

This is the README for fwdpy11_, which is a Python package for forward-time population genetic simulation.  It uses
fwdpp_ as its C++ back-end.

Build status
-----------------------

Master:

.. image:: https://travis-ci.org/molpopgen/fwdpy11.svg?branch=master
    :target: https://travis-ci.org/molpopgen/fwdpy11

Development: 

.. image:: https://travis-ci.org/molpopgen/fwdpy11.svg?branch=dev
    :target: https://travis-ci.org/molpopgen/fwdpy11

Manual
-----------------------

Latest/master:

.. image:: https://readthedocs.org/projects/fwdpy11/badge/?version=latest
	:target: http://fwdpy11.readthedocs.io/en/latest/?badge=latest

Development branch:

.. image:: https://readthedocs.org/projects/fwdpy11/badge/?version=dev
	:target: http://fwdpy11.readthedocs.io/en/dev/?badge=dev

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

fwdpy11 is written for Python 3.  Given time/energy, we *may* look into supporting Python 2.7.  However, it is impossible to support all fwdpy11 features in Python 2.7.  Specifically, it is not possible to provide full pickling support.

.. code-block:: bash

    git submodule init
    git submodule update
    python setup.py build_ext -i
    python test.py

Installation
---------------------------------

Using pip on Linux (or pip3 as appropriate for your system):

.. code-block:: bash

    pip install --upgrade fwdpy11

.. _fwdpy11: http://molpopgen.github.io/fwdpy11
.. _fwdpp: http://molpopgen.github.io/fwdpp
.. _GSL: http://gnu.org/software/gsl
.. _pybind11: https://github.com/pybind/pybind11
.. _multiprocessing: https://docs.python.org/3/library/multiprocessing.html
.. _concurrent.futures: https://docs.python.org/3/library/concurrent.futures.html
