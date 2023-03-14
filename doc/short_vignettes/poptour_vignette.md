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

(poptour_vignette)=

# Properties of populations

This vignette is "getting ahead" of ourselves a bit.
In order for it to be useful, it is best to have a simulated population in hand.
The following hidden code block gives us one.
Feel free to expand it to take a look--the details will make sense once you have studied the later vignettes on setting up and running simulations.
The simulation done here is taken from {ref}`another vignette <gss_vignette>`.

```{code-cell} python
---
tags: ['hide-input']
---
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

p = {
    "nregions": [],
    "gvalue": fwdpy11.Additive(2.0, GSSmo),
    "sregions": [fwdpy11.GaussianS(0, 1., 1, 0.1)],
    "recregions": [fwdpy11.PoissonInterval(0, 1., rho / float(4 * pop.N))],
    "rates": (0.0, 1e-3, None),
    "prune_selected": False,
    "demography": None,
    "simlen": 10 * pop.N,
}
params = fwdpy11.ModelParams(**p)

fwdpy11.evolvets(rng, pop, params, 100, suppress_table_indexing=True)
```

## Basic properties

```{code-cell}
print(f"Total number of diploids = {pop.N}")
print(f"Current generation = {pop.generation}")
```

## Mutations

Mutations are instances of {class}`fwdpy11.Mutation` stored in {attr}`fwdpy11.DiploidPopulation.mutations`.
These objects do not have nice representations in Python:

```{code-cell} python
[i for i in pop.mutations[:5]]
```

So, let's just look at the fields of the first mutation:

```{code-cell}
m = pop.mutations[0]
print(f"position = {m.pos}\neffect_size = {m.s}\ndominance = {m.h}\norigin time = {m.g}")
```

{class}`fwdpy11.Mutation` has other attributes that are not relevant to the type of simulation done here.

## Diploids

Diploids are instances of {class}`fwdpy11.DiploidGenotype`.
They are stored in `fwdpy11.DiploidPopulation.diploids`.
Diploids store the indexes of their individual genomes, which we describe below.

```{code-cell} python
[i for i in pop.diploids[:5]]
```

## Diploid meta data

Diploid individuals have associated data, stored in {attr}`fwdpy11.DiploidPopulation.diploid_metadata`.
The meta data are instances of {class}`fwdpy11.DiploidMetadata`.
That class is also a {class}`numpy.dtype`, allowing us to access the raw data efficiently as a {class}`numpy.recarray`:

```{code-cell} python
import numpy as np
md = np.array(pop.diploid_metadata, copy=False)
md[:5]
```

The field names of the record array exactly match the attribute names of {class}`fwdpy11.DiploidMetadata`.

The record arrays allow efficient calculation of important quantities:

```{code-cell}
print(f"Mean trait value = {md['g'].mean():0.4f}.\nMean fitness = {md['w'].mean():0.4f}")
```

### Ancient samples

When ancient/preserved samples are recorded during simulations, their meta data are stored in :attr:`fwdpy11.DiploidPopulation.ancient_sample_metadata`.

## Haploid genomes

Haploid genomes are instances of {class}`fwdpy11.HaploidGenome`.
These objects contain indexes to mutations.
Let's look at the effect sizes and origin times of mutations in the first non-empty genome:

```{code-cell}
def mut_info(pop, i):
    rv = False
    for m in pop.haploid_genomes[i].smutations:
        print(f"effect size = {pop.mutations[m].s:0.4f}, origin time = {pop.mutations[m].g}")
        rv = True
    return rv
    
for i in pop.diploids:
    if mut_info(pop, i.first) is True:
        break
    elif mut_info(pop, i.second) is True:
        break
```

### Some details

Neutral mutations are never added to haploid genomes.

## Tables

A population contains an instance of {class}`fwdpy11.TableCollection` which is used to represent the genetic ancestry of the sample using the methods described in {cite}`Kelleher2018-mc`.

The table contents are described in the class documentation.
All tables can be accessed either as Python objects or as {class}`numpy.recarray` objects as we say above for diploid meta data.
