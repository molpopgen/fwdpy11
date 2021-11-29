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

(recapitation)=

# Finishing a simulation with msprime

This example is complementary to {ref}`precapitation`.
Rather than starting with a tree sequence from `msprime`, we instead finish a simulation by "coalescing back" the first generation of the simulation using `msprime`.
{cite}`Haller2019-mn` refer to this procedure as "recapitation" of a tree sequence.
In order for recapitation to work correctly, we must pass `preserve_first_generation=True` to {func}`fwdpy11.evolvets`.

First, we'll simulate a population for 10 generations:


```{code-cell} python
import fwdpy11

pop = fwdpy11.DiploidPopulation(100, 1000.0)

pdict = {"rates": (0.0, 0.0, 0.0),
         "gvalue": fwdpy11.Multiplicative(2.),
         "simlen": 10,
        }
params = fwdpy11.ModelParams(**pdict)
rng = fwdpy11.GSLrng(54321)

# We must preserve the founder generation:
fwdpy11.evolvets(rng, pop, params, 100, preserve_first_generation=True)

assert pop.generation == 10

```

Now we use `msprime` to coalesce the founder generation roots back to a common ancestor:

```{code-cell} python
import msprime

# Convert to tskit format
ts = pop.dump_tables_to_tskit()

num_roots_pre_recapitation = ts.first().num_roots

recapitated_ts = msprime.simulate(from_ts=ts)

print(num_roots_pre_recapitation, recapitated_ts.first().num_roots)
```

