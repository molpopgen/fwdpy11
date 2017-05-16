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

fwdpy11 is written for Python 3.  We will not modify the package to be compatibly with Python 2.7.

.. code-block:: bash

    git submodule init
    git submodule update
    python setup.py build_ext -i
    python test.py

.. note::
    The clang compiler is the assumed default on OS X.  However, life is simpler
    if you use gcc.  The setup.py takes a --gcc option that eliminates OS X-specific
    (really Xcode clang-specific) features so that an OS X/gcc build is possible.

Installation
---------------------------------

Using pip on OS X and Linux (or pip3 as appropriate for your system):

.. code-block:: bash

    pip install --upgrade fwdpy11

If you prefer a pip install on OS X using GCC instead of clang:

.. code-block:: bash

    pip install --upgrade fwdpy11 --install-option=--gcc

You may or may not need to prefix the above with

.. code-block:: bash

    CC=gcc CXX=g++

depending on wheter or not your user's `$PATH` is set up to override Xcode's symlink of gcc to clang.

Enabling assertions in the C++ code
------------------------------------------------------------------

The C++ code uses C's assert macros in several places.  These are disabled by default.  However, it can be useful to
enable them when hacking the code.  To do so:

.. code-block:: bash

    python setup.py build_ext -i --debug

.. note::
    Never install the package compiled in debug mode!  First, things will run much more slowly.  
    Second, triggering an assertion will cause the Python interpreter to crash.  These assertions
    exist as a brute-force method to help developers quickly identify bugs.

Bioconda
=================================

fwdpy11 is available through bioconda_ for Linux and for OS X:

.. code-block:: bash

    conda install -c bioconda fwdpy11

The OS X build is built using gcc.

.. _fwdpy11: http://molpopgen.github.io/fwdpy11
.. _fwdpp: http://molpopgen.github.io/fwdpp
.. _GSL: http://gnu.org/software/gsl
.. _pybind11: https://github.com/pybind/pybind11
.. _multiprocessing: https://docs.python.org/3/library/multiprocessing.html
.. _concurrent.futures: https://docs.python.org/3/library/concurrent.futures.html
.. _bioconda: https://bioconda.github.io/
