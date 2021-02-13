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

:::{note}

The details of this section will change once `msprime` 1.0 is released.

:::

### A single deme

```{code-cell} python
import msprime

N = 100
L = 1000.0

ts = msprime.simulate(2*N, Ne = N, length=L, recombination_rate = 1e-3, random_seed = 42)
pop = fwdpy11.DiploidPopulation.create_from_tskit(ts)
```

```{code-cell}
print(f"N = {pop.N}, L = {pop.tables.genome_length}")
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

We can initialize from a multi-deme `msprime` simulation as well:

```{code-cell} python
config = [
    msprime.PopulationConfiguration(sample_size = 50),
    msprime.PopulationConfiguration(sample_size = 150),
]
migration_matrix = [[0.0, 0.1], [0.1, 0.0]]
ts = msprime.simulate(Ne = N,
    population_configurations=config,
    migration_matrix=migration_matrix,
    length = L,
    recombination_rate = 1e-3)
pop = fwdpy11.DiploidPopulation.create_from_tskit(ts)
```

```{code-cell} python
np.unique(np.array(pop.diploid_metadata, copy=False)['deme'],
          return_counts=True)
```

### Limitations and caveats

* Mutations from `msprime` are ignored.
* See {ref}`msprime-subtleties` for more information.
