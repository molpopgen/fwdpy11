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

(mutationdominance_vignette)=

# Mutations with variable dominance

:::{versionadded} 0.13.0
:::

The heterozygous effect of a mutation does not need to be constant.
You may use instances of classes derived from {class}`fwdpy11.MutationDominance` to assign functions to generate the dominance of mutations.

The available classes are:

* {class}`fwdpy11.FixedDominance`
  This class is equivalent to passing in a {class}`float` to the `h` `kwarg`.
  See {ref}`here <des_vignette>`
* {class}`fwdpy11.ExponentialDominance`
* {class}`fwdpy11.UniformDominance`
* {class}`fwdpy11.LargeEffectExponentiallyRecessive`

## Example

```{code-cell}
import fwdpy11

des = fwdpy11.GaussianS(beg=0, end=1, weight=1, sd=0.1,
    h=fwdpy11.LargeEffectExponentiallyRecessive(k=5.0))
```

If we apply this `des` object to a model of a quantitative trait evolving to a sudden "optimum shift", then we see that the larger effect variants present at the end of the simulation to indeed have smaller dominance coefficients.

```{code-cell} python
---
tags: ['hide-input']
---
pop = fwdpy11.DiploidPopulation(500, 1.0)

rng = fwdpy11.GSLrng(54321)

GSSmo = fwdpy11.GSSmo(
    [
        fwdpy11.Optimum(when=0, optimum=0.0, VS=1.0),
        fwdpy11.Optimum(when=10 * pop.N, optimum=1.0, VS=1.0),
    ]
)

rho = 1000.

p = {
    "nregions": [],
    "gvalue": fwdpy11.Additive(2.0, GSSmo),
    "sregions": [des],
    "recregions": [fwdpy11.PoissonInterval(0, 1., rho / float(4 * pop.N))],
    "rates": (0.0, 1e-3, None),
    "prune_selected": False,
    "demography": fwdpy11.DiscreteDemography(),
    "simlen": 10 * pop.N + 200,
}
params = fwdpy11.ModelParams(**p)

fwdpy11.evolvets(rng, pop, params, 100, suppress_table_indexing=True)
```

```{code-cell}
---
tags: ['hide-input']
---
import matplotlib.pyplot as plt
esize = [pop.mutations[m.key].s for m in pop.tables.mutations]
h = [pop.mutations[m.key].h for m in pop.tables.mutations]

f, ax = plt.subplots()
ax.scatter(esize, h)
ax.set_xlabel("Effect size of mutation on trait")
ax.set_ylabel("Dominance")
plt.show();
```

## Using discrete distributions

{class}`fwdpy11.DiscreteDESD` specifies a *D*iscrete *E*ffect *S*ize and *D*ominance joint distribution.
A {class}`list` of {class}`tuple` of `(effect size, dominance, weight)` specify the joint distribution.
For example:

```{code-cell}
import math
import numpy as np

joint_dist = []
for s in np.arange(0.1, 1, 0.1):
    joint_dist.append((-s, math.exp(-s), 1./s))
```

```{code-cell}
des = fwdpy11.DiscreteDESD(beg=0, end=1, weight=1, joint_dist=joint_dist)
print(des.asblack())
```

The result is that mutations with smaller effect sizes are more common (larger weights) and more dominant.
