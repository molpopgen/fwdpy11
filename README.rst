fwdpy11
*************************

This is the README for fwdpy11_, which is a Python package for forward-time population genetic simulation.  It uses
fwdpp_ as its C++ back-end.


.. image:: https://anaconda.org/bioconda/fwdpy11/badges/license.svg
        :target: https://anaconda.org/bioconda/fwdpy11

.. image:: https://anaconda.org/bioconda/fwdpy11/badges/installer/conda.svg
        :target: https://conda.anaconda.org/bioconda

.. image:: https://anaconda.org/bioconda/fwdpy11/badges/version.svg   
	:target: https://anaconda.org/bioconda/fwdpy11

.. image:: https://anaconda.org/bioconda/fwdpy11/badges/platforms.svg   
        :target: https://anaconda.org/bioconda/fwdpy11


Build status
-----------------------

Main:

.. image:: https://github.com/molpopgen/fwdpy11/workflows/Tests/badge.svg?branch=main
    :target: https://github.com/molpopgen/fwdpy11/workflows/Tests/badge.svg?branch=main

.. image:: https://github.com/molpopgen/fwdpy11/workflows/UbuntuStressTest/badge.svg?branch=main
    :target: https://github.com/molpopgen/fwdpy11/workflows/UbuntuStressTest/badge.svg?branch=main

Development: 

.. image:: https://github.com/molpopgen/fwdpy11/workflows/Tests/badge.svg?branch=dev
    :target: https://github.com/molpopgen/fwdpy11/workflows/Tests/badge.svg?branch=dev

.. image:: https://github.com/molpopgen/fwdpy11/workflows/UbuntuStressTest/badge.svg?branch=dev
    :target: https://github.com/molpopgen/fwdpy11/workflows/UbuntuStressTest/badge.svg?branch=dev

Miscellaneous
-----------------------

Python code style:

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black

Features
-----------------------

* Pickle-able population objects
* Parallel computation via multiprocessing_ or concurrent.futures_.
* Custom temporal samplers to analyze populations *during* a simulation may be written in pure Python.
* Flexible interface for simulating models with multiple populations.

Documentation
-----------------------

The manual can be found `here <https://molpopgen.github.io/fwdpy11>`_.

Dependencies
-----------------------

The following must be present on your system:

* GSL_. This is a C library.
  It is available via `conda`.
  fwdpy11 requires version 2.3 or greater.
* cmake_.
  This should be installed by `conda` or your favorite package manager.

fwdpy11_ also uses fwdpp_, which is included as a submodule.

License
-----------------------

GPLv3 or later (See COPYING)

Supported Python version
-------------------------------------------------

fwdpy11 is written for Python 3.
We will not modify the package to be compatible with Python 2.7.


Installation
---------------------------------

Building from the git repository:

.. code-block:: bash

    git submodule init
    git submodule update
    python setup.py build_ext -i
    python -m pytest tests

Using pip on macOS and Linux (or pip3 as appropriate for your system):

.. code-block:: bash

    pip install --upgrade fwdpy11

It is possible that the cmake macros to detect the GSL can fail to detect the correct version.
Issues like this are a basic weakness of cmake.
I've seen this in conda environments, where the macro prefers the system version over the newer version in the environment.
To "fix" this, give it a hint:

.. code-block:: bash

    GSL_ROOT_DIR=/path/to/gsl python3 setup.py build_ext -i

macOS
==================================

On Apple's macOS, we strongly encourage that you use conda with their compiler packages:

.. code-block:: bash

    conda install clang_osx-64 clangxx_osx-64

Installing these packages will mean that you can get away from the relatively old versions of these compilers that ship with Xcode.
However, you do need to add the following flag when building the package:

On macOS versions prior to "Mojave":

.. code-block:: bash

    CONDA_BUILD_SYSROOT=/ python3 setup.py build_ext -i

For later versions, you may omit the environment variable prefix.

Do the same for a `pip install` from the source directory.

Windows
========================================

We have heard positive reports of using fwdpy11 on Windows 10 with the Ubuntu subsystem installed.
For such a system, you may use a Linux conda installer and then install fwdpy11 via bioconda_.

The developers do not have access to this platform, but we are keen to hear of any issues.

Caution
==================================

We use the GitHub "release_" mechanism to make stable versions available.
However, GitHub releases to not include the sub-modules, meaning that the releases themselves cannot be used for installation.
(A related irony is that the Zenodo DOI for the releases are somewhat meaningless.)

To install a specific release:

1. Use pip (see above).
   This is the recommended approach if you do not use conda.
2. Install from bioconda.
   This is the recommended approach.
3. Clone the repo, checkout the release, and update submodules:

.. code-block:: bash

    git clone http://github.com/molpopgen/fwdpy11
    cd fwdpy11
    git submodule init
    git submodule update

The latter method is probably the least appealing.

We have a strict policy of putting releases on PyPi and bioconda_.
If there is a release on PyPi but not on bioconda_, then that is because we identified a bug and pushed a new release before the bioconda_ build happened.
It happens.
That's life.

Bioconda
=================================

fwdpy11 is available through bioconda_ for Linux and for macOS:

.. code-block:: bash

    conda install -c bioconda fwdpy11

.. note::

   Please read the bioconda documentation!  
   The order of channels matters.

.. _fwdpy11: https://github.com/molpopgen/fwdpy11
.. _fwdpp: https://github.com/molpopgen/fwdpp
.. _GSL: http://gnu.org/software/gsl
.. _pybind11: https://github.com/pybind/pybind11
.. _multiprocessing: https://docs.python.org/3/library/multiprocessing.html
.. _concurrent.futures: https://docs.python.org/3/library/concurrent.futures.html
.. _bioconda: https://bioconda.github.io/
.. _release: https://github.com/molpopgen/fwdpy11/releases
.. _cmake: https://cmake.org
