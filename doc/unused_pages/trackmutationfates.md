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

(tracking_mutation_fates)=

# Tracking user-specified new mutations

```{code-cell} python
---
tags: ['hide-input']
---

import fwdpy11
import numpy as np
import msprime
```

```{code-cell} python
import fwdpy11.conditional_models
```


```{code-cell} python
---
tags: ['hide-input']
---

def setup(prune_selected=False):
    # Dropping mutations requires existing
    # ancestry, which we can get either
    # from a burn-in or from msprime.
    initial_ts = msprime.sim_ancestry(
        samples=500,
        population_size=500,
        recombination_rate=1e-1,
        random_seed=43215,
        sequence_length=1.0,
    )

    # Build the pop from msprime output
    pop = fwdpy11.DiploidPopulation.create_from_tskit(initial_ts)

    # Set up basic model parameters
    pdict = {
        "recregions": [fwdpy11.PoissonInterval(0, 1, 1e-1)],
        # Here, 2 means that fitness is multiplicative
        # over 1, 1+hs, 1+2s.
        "gvalue": fwdpy11.Multiplicative(2.0),
        "rates": (0, 0, None),
        "prune_selected": False,
        "simlen": 200,
    }
    params = fwdpy11.ModelParams(**pdict)

    return pop, params
```

## Tracking a mutation for a specified number of generations

```{code-cell} python
ALPHA = -10.0
```


```{code-cell} python
rng = fwdpy11.GSLrng(12345)
pop, params = setup()
```

```{code-cell} python
mutation_data = fwdpy11.conditional_models.NewMutationParameters(
    frequency=fwdpy11.conditional_models.AlleleCount(1),
    data=fwdpy11.NewMutationData(effect_size=ALPHA / 2 / pop.N, dominance=1),
    position=fwdpy11.conditional_models.PositionRange(left=0.49, right=0.51),
)
```


```{code-cell} python
output = fwdpy11.conditional_models.track_added_mutation(
    rng, 
    pop,
    params,
    mutation_data,
    when=3,
    until=7,
)
```

When tracking deleterious variants, it is unlikely that they will be around at the end of the simulation:

```{code-cell} python
try:
    print(output.pop.mutations[output.mutation_index])
except IndexError as _:
    print(f"mutation {output.mutation_index} is no longer in the population!") 
```

### Recording all generations of the mutation's sojourn

So, we've lost all the information about this variant.
That's not so useful.
Let's record all generations of its existence as ancient samples:

```{code-cell} python
rng = fwdpy11.GSLrng(12345)
```

```{code-cell} python
output = fwdpy11.conditional_models.track_added_mutation(
    rng, 
    pop,
    params,
    mutation_data,
    when=3,
    until=7,
    sampling_policy=fwdpy11.conditional_models.AncientSamplePolicy.DURATION,
)
```

Now, our mutation is present in nodes in our tree sequence.
Let's try to print it again:

```{code-cell} python
try:
    print(output.pop.mutations[output.mutation_index])
except IndexError as _:
    output.mutation_index(f"mutation {output.mutation_index} is no longer in the population!") 
```

Let's track this variant's frequency at each time point:

```{code-cell}
for time, nodes, _ in output.pop.sample_timepoints(include_alive=False):
    print(time, len(nodes))
    tree_itr = fwdpy11.TreeIterator(output.pop.tables, nodes)
    for tree in tree_itr:
        for mutation in tree.mutations():
            print(time, tree.leaf_counts(mutation.node), mutation.key)
``` 

