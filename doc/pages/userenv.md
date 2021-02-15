# Setting up user environments for installation

## Standard Unix-like environments

This section assumes an Ubuntu-like Linux.
Users of other distros and Apple's macOS should adjust accordingly.

First, install the GNU Scientific Library, cmake, and a C++ compiler:

```{code-block} bash

sudo apt install g++ cmake libgsl-dev

```

Then, install from PyPi via `pip`:

```{code-block} bash

pip3 install fwdpy11

```

Or, if working from the root of a clone of the repository:

```{code-block} bash

pip3 install -f requirements.txt
pip3 install .

```

The above commands will also work within `virtualenv` for users familiar with that tool.

We refrain from giving detailed instructions for macOS at this time.
This platform has always been a moving target, and major recent changes to Xcode and the switch to ARM chips lead us to anticipate strangeness in the coming months.

## Conda

The basic steps are:

1. Install `conda` for your environment.
   Get a 64-bit installer for Python 3.6 or later.
   For most people, `miniconda` is the way to go.
2. Set up the channel order following instructions [here](<http://bioconda.github.io/user/install.html#set-up-channels>).
   At the time of this writing, the channel order is:

```{code-block} bash

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

```

3. Finally, we can install:

```{code-block} bash

conda install fwdpy11

```

:::{note}

If you have already installed `conda`, you may run into funny issues.
For example, you may already have channels set up in a certain order.
You will want to make sure that the order is correct.

The channels will be in your `~/.condarc` and the contents should look like this:

```
channels:
  - conda-forge
  - bioconda
  - defaults
```

:::

### Working in "environments"

:::{note}

The stability of conda installation into environments is a moving target.
We will try to keep this up to date.

:::

Many users prefer to set up different environments for different projects.
Doing so requires some care.
Currently, `conda-forge`, the channel from which `bioconda` pulls dependencies, builds for Python versions not supported by `bioconda`.
Thus, building an environment like this will often lead to failure:

```{code-block} bash

conda create -n myenv -y

```

The problem here is that installing dependencies that are *not* from `bioconda` is likely to pull in Python 3.9, leaving you unable to install `bioconda` packages.
(Currently, `bioconda` only supports up to Python 3.8.)

Unfortunately, creating a new environment based on Python 3.8 will also lead to failure:

```{code-block} bash

conda create -n myenv python=3.8 -y

```

The problem here is that the Python 3.8 build will be upgraded compared to what came with `miniconda`.
For reasons that we are unable to understand, this upgrade leads to run-time problems in some instances.

Thus, the safest thing to do is to clone your `base`:

```{code-block} bash

conda create -n myenv --clone base

```

Now, `myenv` is an exact copy of `base`, and you can `conda install fwdpy11` successfully.
This solution is safest if you never install anything into `base`, leaving it as it was when you installed `miniconda`.

:::{warning}

When working in environments, never `conda update`!
Updating is likely to do *bad things* to dependencies that you won't notice until you run something.

:::

### Conda compilers

If you will develop or use plugins to `fwdpy11` (see {ref}`here <writingplugins>`), then you will need C++ compilers installed as well as `pybind11` and probably `cmake`.
You will also need these tools if you intend to modify the `fwdpy11` code itself.

On Linux:

```{code-block} bash

conda install cmake pybind11 gxx_linux-64

```

On macOS:

```{code-block} bash

conda install cmake pybind11 clangxx_osx-64

```

(sec_requirements_files)=

## Comments on "requirements" files.

Python uses files with names like `requirements.txt` to specify the package environment that is expected to go along with a certain tool such as this one.
Unfortunately, the contents of such files vary quite a bit across projects.

Here, we take the files `requirements.txt` and `doc/requirements.txt` to specify completely reproducible development and production environments.
Within these files, all dependencies are pinned to specific versions, as are the dependencies of those dependencies.
These files specify the exact versions used for testing the package on Ubuntu environments and for building Docker containers.

The files `requirements.in` and `doc/requirements.in` are looser specifications of these environments.
Here, the dependency pinnings are minimal, trying to keep the build environment compatible with the lowest Python version number that is still supported.
Many users working in extant Python environments may want to use the `.in` versions for installing packages.
Some may also want to remove the installation modifications to various packages in those files.

These pinnings will be updated periodically.

People using `conda` to install `fwdpy11` may consider using the fully-pinned files to install dependencies using `pip`.
