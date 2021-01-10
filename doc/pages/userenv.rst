Setting up user environments for installation
====================================================================================

Standard Unix-like environments
******************************************************************

This section assumes an Ubuntu-like Linux.
Users of other distros and Apple's macOS should adjust accordingly.

First, install the GNU Scientific Library, cmake, and a C++ compiler:

.. code-block:: bash

   sudo apt install g++ cmake libgsl-dev

Then, install from PyPi via `pip`:

.. code-block:: bash

   pip3 install fwdpy11

Or, if working from the root of a clone of the repository:

.. code-block:: bash

   pip3 install -f requirements.txt
   pip3 install .

The above commands will also work within `virtualenv` for users familiar with that tool.

We refrain from giving detailed instructions for macOS at this time.
This platform has always been a moving target, and major recent changes to Xcode and the switch to ARM chips lead us to anticipate strangeness in the coming months.

Conda
*********************************

The basic steps are:

1. Install `conda` for your environment.
   Get a 64-bit installer for Python 3.6 or later.
   For most people, `miniconda` is the way to go.
2. Set up the channel order following instructions `here <http://bioconda.github.io/user/install.html#set-up-channels>`_.
   At the time of this writing, the channel order is:

.. code-block:: bash

   conda config --add channels defaults
   conda config --add channels bioconda
   conda config --add channels conda-forge

3. Finally, we can install:

.. code-block:: bash

   conda install fwdpy11

.. note::

   If you have already installed `conda`, you may run into funny issues.
   For example, you may already have channels set up in a certain order.
   You will want to make sure that the order is correct.

   The channels will be in your `~/.condarc` and the contents should look like this::

      channels:
        - conda-forge
        - bioconda
        - defaults


Conda compilers
++++++++++++++++++++++++++++++++++++++++

If you will develop or use plugins to `fwdpy11` (see :ref:`here <writingplugins>`), then you will need C++ compilers installed as well as `pybind11` and probably `cmake`.
You will also need these tools if you intend to modify the `fwdpy11` code itself.

On Linux:

.. code-block:: bash

   conda install cmake pybind11 gxx_linux-64

On macOS:

.. code-block:: bash

   conda install cmake pybind11 clangxx_osx-64

