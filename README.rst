fwdpy11
*************************

This is the README for fwdpy11_, which is a Python package for forward-time population genetic simulation.  It uses
fwdpp_ as its C++ back-end.


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

Conda status
-----------------------

.. image:: https://anaconda.org/bioconda/fwdpy11/badges/version.svg   
	:target: https://anaconda.org/bioconda/fwdpy11

.. image:: https://anaconda.org/bioconda/fwdpy11/badges/platforms.svg   
        :target: https://anaconda.org/bioconda/fwdpy11

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

License
-----------------------

GPLv3 or later (See COPYING)

Supported Python version
-------------------------------------------------

fwdpy11 is written for Python 3.
We will not modify the package to be compatible with Python 2.7.


Dependencies and installation
---------------------------------

These topics are covered in the user manual:

* People wishing to run `fwpdy11` should see `this section <https://molpopgen.github.io/fwdpy11/pages/userenv.html>`_.
* Those who need to build the package from source should look `the developer's guide <https://molpopgen.github.io/fwdpy11/misc/developersguide.html>`_.


.. _fwdpy11: https://github.com/molpopgen/fwdpy11
.. _fwdpp: https://github.com/molpopgen/fwdpp
.. _GSL: http://gnu.org/software/gsl
.. _pybind11: https://github.com/pybind/pybind11
.. _multiprocessing: https://docs.python.org/3/library/multiprocessing.html
.. _concurrent.futures: https://docs.python.org/3/library/concurrent.futures.html
.. _bioconda: https://bioconda.github.io/
.. _release: https://github.com/molpopgen/fwdpy11/releases
.. _cmake: https://cmake.org
