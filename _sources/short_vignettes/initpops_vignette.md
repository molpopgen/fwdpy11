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

(initpops_vignette)=

# Initializing a population

## Initializing with no previous ancestry

### A single deme

To initialize a {class}`fwdpy11.DiploidPopulation` with $N = 100$ diploids and a genome length of $1000$:

```{code-cell} python
import fwdpy11

pop = fwdpy11.DiploidPopulation(100, 1000.0)
```

We can access these fields:

```{code-cell} python
print(f"N = {pop.N}, L = {pop.tables.genome_length}")
```

The genome length is a property of the {class}`fwdpy11.TableCollection` that records the population's ancestry.

### Multiple initial demes

To initialize $N = 100$ with two demes of sizes 25 and 50, respectively, and the same genome length as above:

```{code-cell} python
pop = fwdpy11.DiploidPopulation([25, 75], 1000.0)
```

The overall $N$ is the same:

```{code-cell} python
print(f"N = {pop.N}")
```

The deme labels are reflected in the individual meta data.
The meta data are a list of instance of {class}`fwdpy11.DiploidMetadata`.
That class is also a {class}`numpy.dtype`, allowing us to access the raw data efficiently as a {class}`numpy.recarray`:

```{code-cell} python
import numpy as np
md = np.array(pop.diploid_metadata, copy=False)
np.unique(md['deme'], return_counts=True)
```

The individuals are intuitively ordered:

```{code-cell} python
print(np.where(md['deme'] == 0)[0].min(),
      np.where(md['deme'] == 0)[0].max(),
      np.where(md['deme'] == 1)[0].min(),
      np.where(md['deme'] == 1)[0].max())
```

(precapitation)=

## Initializing with ancestry from msprime

### A single deme

```{code-cell} python
import msprime

N = 100
L = 1000.0

ts = msprime.sim_ancestry(N, population_size = N, sequence_length=L, recombination_rate = 1e-3, random_seed = 42)
pop = fwdpy11.DiploidPopulation.create_from_tskit(ts)
assert pop.N == 100
assert pop.tables.genome_length == L
```

Now, our tables have some data in them.
Edges are stored as a list of {class}`fwdpy11.Edge` and nodes as a list of {class}`fwdpy11.Node`.
As with meta data, these are {class}`numpy.dtype`s:


```{code-cell} python
np.array(pop.tables.edges, copy=False)[:10]
```

```{code-cell} python
np.array(pop.tables.nodes, copy=False)[:10]
```

### Multiple demes

We can initialize from a multi-deme `msprime` simulation as well.
The simplest way to do this is to use [demes](https://popsim-consortium.github.io/demes-spec-docs/main/tutorial.html) to specify the model:

```{code-cell} python
import demes

model_yaml = """
description:
  Example from the fwdpy11 manual
time_units: generations
demes:
  - name: deme0
    epochs:
      - start_size: 50
        end_time: 0
  - name: deme1
    epochs:
      - start_size: 150
        end_time: 0
migrations:
  - demes: [deme0, deme1]
    rate: 0.1    
"""

graph = demes.loads(model_yaml)

demography = msprime.Demography.from_demes(graph)
ts = msprime.sim_ancestry(
    samples={0: 50, 1: 150},
    demography=demography,
    sequence_length=100,
    recombination_rate=1e-3)
pop = fwdpy11.DiploidPopulation.create_from_tskit(ts)
assert pop.N == 200
ds = pop.deme_sizes(as_dict=True)
assert ds[0] == 50
assert ds[1] == 150
```

### Importing msprime tree sequences with mutations.

* See {ref}`import_mutations_from_tskit_vignette` for more information.
