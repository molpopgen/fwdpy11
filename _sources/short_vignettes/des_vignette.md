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

(des_vignette)=

# Distributions of effect sizes

:::{note}
The lists of objects described here are passed to the `sregions` parameter when initializing instances of {class}`fwdpy11.ModelParams`.
:::

One of the main reasons to perform forward simulations is to be able to model mutations affecting individual fitness.
To do so, we need to specify both mutation rates and the resulting effect sizes.

`fwdpy11` works by specifying an overall mutation rate to variants affecting fitness (see {ref}`here <model-params>`).
Given that a mutation occurs, we need to specify its "effect size".

`fwdpy11` chooses the effect size of a new mutation by first determining what *region* is mutated and then generating a mutation from the distribution of effect size associated with that region.

Each region is represented by instances of classes derived from the ABC {class}`fwdpy11.Sregion`.
Each instance is associated with a *weight*.
These weights establish the relative probability that a mutation comes from a given region.
Thus, given an overall mutation rate to non-neutral variants, instances of "`sregions`" are used to set up a multinomial distribution for generating new mutations.

The following sets up a model where mutations have a constant effect size ({math}`s=-0.01`), dominance {math}`h=0.25`, and occur uniformly on the interval {math}`[0, 1)`: 

```{code-cell} python
import fwdpy11

sregions = [fwdpy11.ConstantS(beg=0.0, end=1.0, weight=1.0, s=-0.01, h=0.25)]
```

The previous example uses argument names for clarity, and the following is equivalent, with the `int` values getting converted to `float` automatically:

```{code-cell} python
sregions = [fwdpy11.ConstantS(0, 1, 1, -0.01, 0.25)]
sregions[0]
```

Note that the constructor parameters for these classes often have default values--see the specific class documentation for details. 

In some scenarios, it is useful to think about the distribution of effect sizes as scaled with respect to the population size.
For example, selection coefficients may be exponentially-distributed with a mean of {math}`2Ns`. 
To do this in `fwdpy11`:

```{code-cell} python
# ALPHA = 2Ns
MEAN_ALPHA = -10
N = 1000
sregions = [fwdpy11.ExpS(0, 1, 1, MEAN_ALPHA, scaling=2 * N)]
sregions[0]
```

## Region weighting

When multiple "sregion" objects are used, the default behavior is to multiply the input `weight` by `end-beg`:

```{code-cell} python
sregions = [
   fwdpy11.ExpS(beg=0.0, end=1.0, weight=1.0, mean=-0.2),
   fwdpy11.ConstantS(beg=1.0, end=3.0, weight=1.0, s=-0.1),
]
sregions
```

Here, the input `weight` is interpreted to mean the weight "per site" is constant.
In this example, twice as many mutations will have positions in {math}`[1, 3)` as from {math}`[0, 1)`. 
To change the default behavior, one can prevent the coupling between input `weight` and region length:

```{code-cell} python
sregions = [
   fwdpy11.ExpS(beg=0.0, end=1.0, weight=1.0, mean=-0.2, coupled=False),
   fwdpy11.ConstantS(beg=1.0, end=3.0, weight=1.0, s=-0.1, coupled=False),
]
sregions
```

The absolute values of the `weight` parameters themselves is irrelevant.
The only thing that matters is the *relative* values from region to region.
Simulations based on the above examples would give the same results if the `weight` were 42 or 73.06.
Therefore, we can recreate our first example with code like the following:

```{code-cell} python
sregions = [
   fwdpy11.ExpS(beg=0.0, end=1.0, weight=56.0, mean=-0.2, coupled=False),
   fwdpy11.ConstantS(beg=1.0, end=3.0, weight=112.0, s=-0.1, coupled=False),
]
sregions
```

In the above example, twice as many mutations occur in the second region because the weights have relative values of 2:1.

:::{note}

Different regions are allowed to overlap, allowing the simulation of concepts like "coding regions" where the DFE are a weighted mixture from multiple distributions, etc.

:::

## Setting the dominance of new mutations

The dominance of a new mutation is set by the `h` parameter during initialization:

```{code-cell} python
fwdpy11.ExpS(beg=0.0, end=1.0, weight=1.0, mean=-0.2, h=1.0)
fwdpy11.ExpS(beg=0.0, end=1.0, weight=1.0, mean=-0.2, h=0.0)
```

## Built-in distributions of effect sizes

* {class}`fwdpy11.ConstantS`
* {class}`fwdpy11.UniformS`
* {class}`fwdpy11.ExpS`
* {class}`fwdpy11.GammaS`
* {class}`fwdpy11.GaussianS`
* {class}`fwdpy11.MultivariateGaussianEffects`
* {class}`fwdpy11.LogNormalS`
* {class}`fwdpy11.DiscreteDESD`

## Labelling mutations from different regions

It may be of use to know what region a mutation came from.
To do so, give a nonzero value to the `label` argument:

```{code-cell} python
fwdpy11.ConstantS(beg=0.0, end=1.0, weight=1.0, s=1e-3, label=1)
```

At the end of the simulation, mutations from this region will have the `label` value stored in the attribute {attr}`fwdpy11.Mutation.label`.

The value of `label` must fit into a 16-bit unsigned integer, *e.g.*, {attr}`numpy.uint16`.
Larger values, or negative values, will result in exceptions.
The following example tries to use a value one larger than the maximum allowed:

```{code-cell} python
---
tags: [raises-exception]
---
import numpy as np

fwdpy11.ConstantS(
    beg=0.0, end=1.0, weight=1.0, s=1e-3, label=np.iinfo(np.uint16).max + 1
)
```

