fwdpy11
*************************

This is the README for fwdpy11_, which is a Python package for forward-time population genetic simulation.  It uses
fwdpp_ as its C++ back-end.

Dependencies
-----------------------

The following must be present on your system:

* GSL_
* pybind11_

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
