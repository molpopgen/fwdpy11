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
         "demography": fwdpy11.ForwardDemesGraph.tubes([100], burnin=10, burnin_is_exact=True),
        }
params = fwdpy11.ModelParams(**pdict)
rng = fwdpy11.GSLrng(54321)

# We must preserve the founder generation:
fwdpy11.evolvets(rng, pop, params, 100, preserve_first_generation=True)

assert pop.generation == 10

```

Now we use `msprime` to coalesce the founder generation roots back to a common ancestor:

```{code-cell} python
import demes
import msprime


# Convert to tskit format
ts = pop.dump_tables_to_tskit()

num_roots_pre_recapitation = ts.first().num_roots

yaml=f"""
time_units: generations
demes:
 - name: deme0
   epochs:
    - start_size: {pop.N}
"""
graph = demes.loads(yaml)
demography = msprime.Demography.from_demes(graph)

recapitated_ts = msprime.sim_ancestry(demography=demography, initial_state=ts)

print(num_roots_pre_recapitation, recapitated_ts.first().num_roots)
```
## Important considerations

The previous example was very simple.
The model involved no recombination and a single deme.

### Demography

In the above example, we had to provide `msprime` a demographic model.
That model must be the correct model for the first generation of your simulation!

### Recombination/genetic maps

You must take care to proved `msprime` with a correct genetic map!
This software and `msprime` differ in some key ways:

* Here, rates are per *genomic segment*.  In `msprime`, rates are per "base pair".
* For forward simulations with unlinked regions, you must take special care when defining a recombination map in `msprime`

### Using the proper backwards-time model

`msprime` supports a few different models of the backwards process.
The two most relevant to this discussion are the "Hudson" and "discrete-time Wright-Fisher" models.
The former is the continuous approximation with recombination and is what most people think of when they think "coalescent simulation".
The latter model allows you to couple Wright-Fisher dynamics to model the recent past with the Hudson algorithm to simulate more ancient events.
You will want to read the msprime [documentation](https://tskit.dev/msprime/docs/stable/intro.html) and the literature cited therein to make a decision about how to best model the ancestry of the ancestral roots of your simulation.


