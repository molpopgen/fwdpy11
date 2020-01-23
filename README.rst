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

Master:

.. image:: https://travis-ci.org/molpopgen/fwdpy11.svg?branch=master
    :target: https://travis-ci.org/molpopgen/fwdpy11

.. image:: https://circleci.com/gh/molpopgen/fwdpy11/tree/master.svg?style=svg
    :target: https://circleci.com/gh/molpopgen/fwdpy11/tree/master

Development: 

.. image:: https://travis-ci.org/molpopgen/fwdpy11.svg?branch=dev
    :target: https://travis-ci.org/molpopgen/fwdpy11

.. image:: https://circleci.com/gh/molpopgen/fwdpy11/tree/dev.svg?style=svg
    :target: https://circleci.com/gh/molpopgen/fwdpy11/tree/dev

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

* GSL_. This is a C library.  It is available via `conda`.  fwdpy11 requires version 2.3 or greater.
* pybind11_. This should be installed `conda` as appropriate for your system, or via your system's package manager or
  manually.  See note below.
* cmake_. This should be installed by `conda` or your favorite package manager.

fwdpy11_ also uses fwdpp_, which is included as a submodule.

.. note::

    The C++ modules are built using cmake_, which requires that pybind11_'s cmake macros are visible.
    Installing pybind11_ via `pip` does **not** install the macros.  However, installs using `conda`, 
    apt-get, or manual installation from source will install both the Python module and the cmake macros.

License
-----------------------

GPLv3 or later (See COPYING)

Suppored Python version
-----------------------

fwdpy11 is written for Python 3.  We will not modify the package to be compatible with Python 2.7.


Installation
---------------------------------

Building from the git repository:

.. code-block:: bash

    git submodule init
    git submodule update
    python setup.py build_ext -i
    python -m unittest discover tests

Using pip on macOS and Linux (or pip3 as appropriate for your system):

.. code-block:: bash

    pip install --upgrade fwdpy11

It is possible that the cmake macros to detect the GSL can fail to detect the correct version.  Issues like this are a
basic weakness of cmake.  I've seen this in conda environments, where the macro prefers the system version over the
newer version in the environment.  To "fix" this, give it a hint:

.. code-block:: bash

    GSL_ROOT_DIR=/path/to/gsl python3 setup.py build_ext -i

macOS
==================================

On Apple's macOS, we strongly encourage that you use conda with their compiler packages:

.. code-block:: bash

    conda install clang_osx-64 clangxx_osx-64

Installing these packages will mean that you can get away from the relatively old versions of these compilers that ship
with Xcode.   However, you do need to add the following flag when building the package:

On macOS versions prior to "Mojave":

.. code-block:: bash

    CONDA_BUILD_SYSROOT=/ python3 setup.py build_ext -i

For later versions, you may omit the environment variable prefix.

Do the same for a `pip install` from the source directory.

Windows
========================================

We have heard positive reports of using fwdpy11 on Windows 10 with the Ubuntu subsystem installed.  For such
a system, you may use a Linux conda installer and then install fwdpy11 via bioconda_.

The developers do not have access to this platform, but we are keen to hear of any issues.

Caution
==================================

We use the GitHub "release_" mechanism to make stable versions available.  However, GitHub releases to not include the
sub-modules, meaning that the releases themselves cannot be used for installation.  (A related irony is that the Zenodo
DOI for the releases are somewhat meaningless.)

To install a specific release:

1. Use pip (see above).  This is the recommended approach if you do not use conda.
2. Install from bioconda.  This is the recommended approach.
3. Clone the repo, checkout the release, and update submodules:

.. code-block:: bash

    git clone http://github.com/molpopgen/fwdpy11
    cd fwdpy11
    git submodule init
    git submodule update

The latter method is probably the least appealing.

We have a strict policy of putting releases on PyPi and bioconda_.  If there is a release on PyPi but not on bioconda_,
then that is because we identified a bug and pushed a new release before the bioconda_ build happend.  It happens.
That's life.

Enabling code profiling
-------------------------------------------------------------------

By default, fwdpy11 is compiled with aggressive optimizations to help reduce the library size. One side effect
is that it becomes impossible to accurately profile the code.  To override these defaults:

.. code-block:: bash

   python setup.py build_ext -i --enable-profiling

.. note::

   The package should not be installed with profiling enabled. This method of building
   is for developers who need to accurately profile the C++ back-end.  Also note that
   only the main package is affected.  Building the unit test modules is not affected.

Disabling link-time optimization (LTO)
------------------------------------------------------------------

LTO is enabled by default and reduced the final library size substantially. However, it takes a
long time and is therefore a drag during development.  To disable it:

.. code-block:: bash

   python setup.py build_ext -i --disable_lto

.. note::

   This option only affects the main package and not the unit tests.


Enabling debugging symbols in the C++ code
------------------------------------------------------------------

.. code-block:: bash

    python setup.py build_ext -i --debug

Debug mode disables all compiler optimizations, allows C-like assertions, and generated debug symbols.

.. note::
    Never install the package compiled in debug mode!  First, things will run much more slowly.  
    Second, triggering an assertion will cause the Python interpreter to crash.  These assertions
    exist as a brute-force method to help developers quickly identify bugs.

Enabling assertions in the C++ code
------------------------------------------------------------------

The fwdpp library code uses C's assert macros in several places.  These are disabled by default.  However, it can be useful to
enable them when hacking the code.  To do so, you must manually set your compiler flags with cmake:

.. code-block:: bash
    
    cmake . -DCMAKE_CXX_FLAGS="-UNDEBUG -O2 -g"

When compiling this way, fwdpy11 makes some extra checks that will throw `RuntimeError` if they fail.  The fwdpp_ back
end also makes extra checks.  If those fail, `abort` will be called, which will crash the Python interpreter.  Thus,
compiling with this option is a "serious debugging mode only" option.

Enabling aggressive debugging of C++ STL templates using GCC
------------------------------------------------------------------------------------------------------------------------------------

Use the following flags to enable an "extreme" debugging mode of the C++ standard template library:

.. code-block:: bash

   CXXFLAGS="-D_GLIBCXX_CONCEPT_CHECKS -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC" \
      CPPFLAGS="-D_GLIBCXX_CONCEPT_CHECKS -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC" python3 setup.py build_ext -i

Bioconda
=================================

fwdpy11 is available through bioconda_ for Linux and for macOS:

.. code-block:: bash

    conda install -c bioconda fwdpy11

.. _fwdpy11: http://molpopgen.github.io/fwdpy11
.. _fwdpp: http://molpopgen.github.io/fwdpp
.. _GSL: http://gnu.org/software/gsl
.. _pybind11: https://github.com/pybind/pybind11
.. _multiprocessing: https://docs.python.org/3/library/multiprocessing.html
.. _concurrent.futures: https://docs.python.org/3/library/concurrent.futures.html
.. _bioconda: https://bioconda.github.io/
.. _release: https://github.com/molpopgen/fwdpy11/releases
.. _cmake: https://cmake.org
