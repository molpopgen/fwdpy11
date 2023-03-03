---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

(timeseries)=

# Time series analysis

fwdpy11 allows populations to be analyzed during the course of a simulation.  Any valid callable type
may be used.  The most useful such type will be a class defining `__call__` to update member data.

There are built-in types to handle common scenarios:

```{eval-rst}
.. autoclass:: fwdpy11.RandomAncientSamples
    :members:

    This type randomly samples diploids at predetermined time points.
    For example, if we wanted to randomly sample 10 diploids
    every N generations for the first 9N generations of a simulation:
```

```{code-cell} python
import fwdpy11
import numpy as np

N = 1000
seed = 42
samplesize = 10
r = fwdpy11.RandomAncientSamples(seed, samplesize, np.arange(0, 10*N, N))
```

When implementing your own callables, the arguments to the type must be of the following form:

```{code-block} python

def timeseries_fxn(pop, sampler):
    pass

```

The first argument is the population itself.  The second is an instance of the following class:

```{eval-rst}
.. autoclass:: fwdpy11.SampleRecorder

    Instances of this type are *never* created by a user.  Rather,
    they are generated internally and passed to user-defined
    callable types.  These types are responsible for passing
    in the indexes of individuals to be preserved as
    samples in the tree sequences.  Two member functions exist
    for such recording.  The first function adds samples one
    at a time:

    .. autofunction:: fwdpy11.SampleRecorder.add_sample

    The second is more efficient, adding in "batch" mode
    via a numpy array with dtype uint32:

    .. autofunction:: fwdpy11.SampleRecorder.assign

    The following is essentially pseudocode to illustrate the process:
```

```{code-cell} python
import fwdpy11
s = fwdpy11.SampleRecorder()
s.add_sample(10)
s.assign(np.arange(10, dtype=np.uint32))
```

The important thing to keep in mind is that user-defined callable
types can do essentially *anything*.  
When recording samples
to preserve, "anything goes".  For example, you could pick samples with the top `x` and
the lowest `y` fitnesses.  You could pick individuals based on their genetic load, whether or
not they carry a specific variant, etc..  It is up to you.
