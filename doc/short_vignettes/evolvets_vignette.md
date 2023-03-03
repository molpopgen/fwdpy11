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

(evolvets_vignette)=

# Evolving with tree sequence recording

:::{note}
The function described in this section has many `kwargs` specifying various advanced features.
Some of these will be the subject of later vignettes.

Future releases are likely to change these arguments to be `kwarg` only.
It is also possible that we will move them into a class that can be pickled, added to `tskit` tree files, etc..
:::

The function {func}`fwdpy11.evolvets` evolves a population with tree sequence recording.

## Our model

By way of example, let's use the following model:

```{code-cell}
import fwdpy11

pop = fwdpy11.DiploidPopulation(500, 1.0)

rng = fwdpy11.GSLrng(54321)

GSSmo = fwdpy11.GSSmo(
    [
        fwdpy11.Optimum(when=0, optimum=0.0, VS=1.0),
        fwdpy11.Optimum(when=10 * pop.N - 200, optimum=1.0, VS=1.0),
    ]
)

rho = 1000.

des = fwdpy11.GaussianS(beg=0, end=1, weight=1, sd=0.1,
    h=fwdpy11.LargeEffectExponentiallyRecessive(k=5.0))

p = {
    "nregions": [],
    "gvalue": fwdpy11.Additive(2.0, GSSmo),
    "sregions": [des],
    "recregions": [fwdpy11.PoissonInterval(0, 1., rho / float(4 * pop.N))],
    "rates": (0.0, 1e-3, None),
    "prune_selected": False,
    "demography": None,
    "simlen": 10 * pop.N,
}
params = fwdpy11.ModelParams(**p)
```

The tables recording the genetic ancestry of the population are empty.
For example, there are no edges:

```{code-cell}
len(pop.tables.edges)
```

Now, we can call our function to evolve our population.
We will apply the table simplification algorithm of {cite}`Kelleher2018-mc` every 100 generations.

```{code-cell}
fwdpy11.evolvets(rng, pop, params, simplification_interval=100) 
```

Now, the edge table is populated:

```{code-cell}
import numpy as np

np.array(pop.tables.edges, copy=False)[:10]
```


