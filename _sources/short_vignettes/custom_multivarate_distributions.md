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

(custom_multivariate_distributions_vignette)=

# "Custom" multivariate distributions

```{code-cell} python
---
tags: ['hide-input']
---
import demes
import demesdraw
import fwdpy11
import numpy as np
```

The previous two sections cover cases where the methods for generating
deviates from a multivariate distribution are straightforward and agreed
upon.

In order to simulate multivariate distributions of effect sizes based on
{class}`fwdpy11.Sregion` types, we follow a fairly intuitive approach
described in {cite}`Xue-Kun_Song2000-qn`.  Briefly, the multivariate Gaussian kernel is
used to produce deviates.  Then, the quantiles from the cummulative distribution
of each marginal Gaussian are used to generate a deviate from the desired output distribution of interest.

For a simulation with `n` populations we need:

* A {class}`list` of `n` {class}`fwdpy11.Sregion` objects
* An array of `n` means for the multivariate Gaussian
* An `n-by-n` covariance matrix for the multivariate
  Gaussian

We will use the same demographic model as in the previous vignettes.
The details can be viewed by expanding the "click to show" button below".

```{code-cell}python
---
tags: ['hide-input']
---
yaml = """
description: Island model forever
time_units: generations
demes:
  - name: A
    epochs:
     - start_size: 100
  - name: B
    epochs:
     - start_size: 100
migrations:
  - demes: [A, B]
    rate: 0.10
"""
g = demes.loads(yaml)
model = fwdpy11.discrete_demography.from_demes(g, burnin=1)
initial_sizes = [v for v in model.metadata["initial_sizes"].values()]
pdict = {
    "nregions": [],
    "recregions": [],
    "sregions": None, # Will get filled in below
    "rates": (0, 5e-3, None),
    "demography": model,
    "simlen": model.metadata["total_simulation_length"],
    "gvalue": fwdpy11.Multiplicative(ndemes=2, scaling=2),
}
rng = fwdpy11.GSLrng(123512)
```

The following generates exponentially distributed effect sizes in each deme
with a high correlation across demes:

```{code-cell} python
mvdes = fwdpy11.mvDES(
    [fwdpy11.ExpS(0, 1, 1, -0.5)] * 2,
    np.zeros(2),
    np.matrix([1, 0.9, 0.9, 1]).reshape((2, 2)),
)
pdict["sregions"] = [mvdes]
params = fwdpy11.ModelParams(**pdict)
pop = fwdpy11.DiploidPopulation(initial_sizes, 1.0)
fwdpy11.evolvets(rng, pop, params, 10)
for i in pop.tables.mutations:
    print(pop.mutations[i.key].esizes)
```

We can mix and match our distributions.  Here, the distribution of effect
sizes in deme 0 is exponential and the distribution in deme 1 is gamma.  The
two distributions have means with opposite signs and the magnitudes of the
marginal deviates negatively covary:

```{code-cell} python
mvdes = fwdpy11.mvDES(
    [fwdpy11.ExpS(0, 1, 1, -0.5), fwdpy11.GammaS(0, 1, 1, mean=0.1, shape_parameter=1)],
    np.zeros(2),
    np.matrix([1, -0.9, -0.9, 1]).reshape((2, 2)),
)
pdict["sregions"] = [mvdes]
params = fwdpy11.ModelParams(**pdict)
pop = fwdpy11.DiploidPopulation(initial_sizes, 1.0)
fwdpy11.evolvets(rng, pop, params, 10)
for i in pop.tables.mutations:
    print(pop.mutations[i.key].esizes)
```

The type {class}`fwdpy11.ConstantS` has intuitive behavior:

```{code-cell} python
mvdes = fwdpy11.mvDES(
    [fwdpy11.ExpS(0, 1, 1, -0.5), fwdpy11.ConstantS(0, 1, 1, -0.1)],
    np.zeros(2),
    np.matrix([1, -0.9, -0.9, 1]).reshape((2, 2)),
)
pdict["sregions"] = [mvdes]
params = fwdpy11.ModelParams(**pdict)
pop = fwdpy11.DiploidPopulation(initial_sizes, 1.0)
rng = fwdpy11.GSLrng(1010)
fwdpy11.evolvets(rng, pop, params, 10)
for i in pop.tables.mutations:
    print(pop.mutations[i.key].esizes)
```
