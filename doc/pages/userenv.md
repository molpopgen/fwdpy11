# Setting up user environments for installation

## Linux and macOS

For those using 64-bit Linux systems and macOS on Intel chips, `fwdpy11` can be installed directly from [PyPi](https://pypy.org) with no need to compile anything.
The recommended procedure is to install into a virtual environment:

```{code-block} bash
mkdir ~/venvs
python3 -m venv ~/venvs/fwdpy11
source activate ~/venvs/fwdpy11/bin/activate
python -m pip install --upgrade pip
python -m pip install fwpdy11
```

:::{note}

The reason to upgrade `pip` is so that dependencies are properly-handled.
This step will not be required for all users, but we recommend it just in case

:::

If you must install from source, see the {ref}`developer's guide <developersguide>`.

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
