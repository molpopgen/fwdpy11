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

(import_mutations_from_tskit_vignette)=

# Importing mutations from tskit

```{code-cell} python
import msprime

ts = msprime.sim_ancestry(1000,
                          recombination_rate=1e-9,
                          sequence_length=1e6,
                          population_size=1000,
                          random_seed=54321)
```

Make a copy of the tables to work with:

```{code-cell} python
tables = ts.tables
```

Add the `fwdpy11` mutation metadata schema:

```{code-cell} python
import fwdpy11

tables.mutations.metadata_schema = fwdpy11.tskit_tools.metadata_schema.MutationMetadata
```

Get a list of candidate edges:

```{code-cell} python
edge_indexes = []
edge_weights = []

for i in range(0, tables.edges.num_rows):
    edge = tables.edges[i]
    ptime = tables.nodes.time[edge.parent]
    ctime = tables.nodes.time[edge.child]
    # We must assign mutations integer-valued
    # times to be compatible with discrete-time
    # forward simulations.  To do so, we filter
    # on branches greater than a generation long.
    if ptime - ctime > 1.0:
        edge_indexes.append(i)
        edge_weights.append((ptime-ctime)*(edge.right - edge.left))
```


We will add 10 mutations to our table collection:

* Effect sizes will be `N(0, 0.1)`.
* Dominance will be 1.

First, choose 10 random edges:

```{code-cell} python
import numpy as np
np.random.seed(666 * 42)

edge_weights = np.array(edge_weights)
edge_weights /= np.sum(edge_weights)

chosen_edges = np.random.choice(edge_indexes, size=10, p=edge_weights)
```

::: {warning}
The method used here to obtain candidate branches for placing mutations and sampling them is **absolutely the wrong thing to do**.
The code shown above is for illustration purposes here.
The approach used here violates many core evolutionary principles!
:::

```{code-cell} python
# We must make sure that mutations
# occur at unique sites because fwdpy11
# is limited to the infinitely-many sites model.
sites = set()
for i in chosen_edges:
    edge = tables.edges[i]

    position = np.random.uniform(edge.left, edge.right, 1)[0]
    while position in sites:
        position = np.random.uniform(edge.left, edge.right, 1)[0]
    sites.add(position)

    site = tables.sites.add_row(position, '0')
    effect_size = np.random.normal(loc=0.0, scale=0.1, size=1)[0]

    # Here, we cheat and take the largest integer value >=
    # the node child time.  In production, something else
    # should be done
    time = np.ceil(tables.nodes.time[edge.child])
    assert time < tables.nodes.time[edge.parent]
    print(time, position, effect_size)

    # Build the mutation metadata
    md = {'s': effect_size,
          'h': 1.0,
          'origin': int(time),
          # NOTE: always write the
          # next 3 lines as shown here.
          # The fwdpy11 back end will do
          # the right thing.
          # A future release will provide a
          # nicer API so that you only need
          # to provide the previous 3 fields.
          'neutral': 0,
          'label': np.uint16(0),
          'key': np.uint64(0)
         }
    tables.mutations.add_row(site, edge.child,
                             '1', time=time,
                             metadata=md)

tables.sort()
ts_with_muts = tables.tree_sequence()
```

Create a population and lift over the mutations:

```{code-cell} python
pop = fwdpy11.DiploidPopulation.create_from_tskit(ts_with_muts, import_mutations=True)
```

Print the number of times each mutation appears in the population and its internal data:

```{code-cell} python
for c, m in zip(pop.mcounts, pop.mutations):
    print(c, "->", m)
```

Set up a model of additive effects on a phenotype and Gaussian stabilizing
selection with an optimum of 1 and the (inverse) strength of stabilizing
selection is also 1.0.

::: {note}
Pay attention to the recombination rate here: it is per genome in `fwdpy11`
but per "base pair" in `msprime`!
:::

```{code-cell} python
pdict = {'nregions': [], 'sregions': [],
         'recregions': [fwdpy11.PoissonInterval(0, pop.tables.genome_length,
                                                1e-9*pop.tables.genome_length)],
         'rates': (0., 0., None),
         'gvalue': fwdpy11.Additive(scaling=2, 
                                     gvalue_to_fitness=fwdpy11.GSS(optimum=1.0, VS=1.0)),
         'simlen': 100,
         'demography': fwdpy11.DiscreteDemography(),
         'prune_selected': False
         }
params = fwdpy11.ModelParams(**pdict)
```

Set up recording some data:

```{code-cell} python
from dataclasses import dataclass

@dataclass
class SimData:
    generation: int
    mean_phenotype: float
    mean_fitness: float

@dataclass
class Recorder:
    data: list

    def __call__(self, pop, _):
        md = np.array(pop.diploid_metadata)
        mean_pheno = md['g'].mean()
        mean_fitness = md['w'].mean()
        self.data.append(SimData(pop.generation, mean_pheno, mean_fitness))
```

Run the simulation:

```{code-cell} python
rng = fwdpy11.GSLrng(2351235)
recorder = Recorder(data = [])
fwdpy11.evolvets(rng, pop, params, simplification_interval=50,
                 recorder=recorder, suppress_table_indexing=True)
assert pop.generation == 100
```

Plot our results:

```{code-cell}
import matplotlib.pyplot as plt

f, ax = plt.subplots()
g = [i.generation for i in recorder.data]
p = [i.mean_phenotype for i in recorder.data]
w = [i.mean_fitness for i in recorder.data]
ax.plot(g, p, label=r'$\bar{g}$')
ax.plot(g, w, label=r'$\bar{w}$')
ax.set_xlabel('Generation')
ax.set_ylabel('Value')
ax.legend(loc='best');
```

What mutations are left at the end of the sim?

```{code-cell} python
for c, m in zip(pop.mcounts, pop.mutations):
    print(c, "->", m)
```
