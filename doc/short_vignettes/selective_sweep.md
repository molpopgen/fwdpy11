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

(selective_sweeps)=

# Selective sweeps

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
import fwdpy11.tskit_tools
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
        recombination_rate=1e-4,
        random_seed=43215,
        sequence_length=10000.0,
    )

    # Build the pop from msprime output
    pop = fwdpy11.DiploidPopulation.create_from_tskit(initial_ts)

    # Set up basic model parameters
    pdict = {
        "recregions": [fwdpy11.PoissonInterval(0, 10000, 1e-4, discrete=True)],
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

## From a new mutation

```{code-cell} python
ALPHA = 1000.0
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
output = fwdpy11.conditional_models.selective_sweep(
    rng, 
    pop,
    params,
    mutation_data,
    fwdpy11.conditional_models.GlobalFixation()
)
```

```{code-cell} python
assert output.pop.generation == params.simlen
assert pop.generation == 0
```

```{code-cell} python
print(output.pop.mutations[output.mutation_index])
```

```{code-cell}
for fixation, time in zip(output.pop.fixations, output.pop.fixation_times):
    print(fixation, time)
```

```{code-cell}
FIXATION_TIME = output.pop.fixation_times[0]
```



### Recording the generation when fixation happened

```{code-cell} python
rng = fwdpy11.GSLrng(12345)
```

```{code-cell} python
output = fwdpy11.conditional_models.selective_sweep(
    rng, 
    pop,
    params,
    mutation_data,
    fwdpy11.conditional_models.GlobalFixation(),
    sampling_policy=fwdpy11.conditional_models.AncientSamplePolicy.COMPLETION,
)
```

```{code-cell} python
assert len(output.pop.ancient_sample_nodes) == 2 * output.pop.N
assert output.pop.fixation_times[output.mutation_index] == FIXATION_TIME
```


```{code-cell} python
node_array = np.array(output.pop.tables.nodes, copy=False)
ancient_sample_node_times = \
    node_array["time"][output.pop.ancient_sample_nodes]
assert np.all([ancient_sample_node_times == \
    output.pop.fixation_times[output.mutation_index]])
```

## From a standing variant

The recipes for a standing variant are identical to those show above, except that one uses {class}`fwdpy11.conditional_models.AlleleCountRange` or {class}`fwdpy11.conditional_models.FrequencyRange` to specify the starting frequencies.
