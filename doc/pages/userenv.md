# Setting up user environments for installation

## Linux and macOS

For those using 64-bit Linux systems and macOS on Intel chips, `fwdpy11` can be installed directly from [PyPi](https://pypy.org) with no need to compile anything.
The recommended procedure is to install into a virtual environment:

```{code-block} bash
mkdir ~/venvs
python3 -m venv ~/venvs/fwdpy11
source activate ~/venvs/fwdpy11/bin/activate
python -m pip install --upgrade pip
pip install fwpdy11
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
