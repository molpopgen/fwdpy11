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

(tskitconvert_vignette)=

# Converting data to tskit format

```{code-cell} python
import msprime
import fwdpy11

# Random number generators are initialized with seeds
rng = fwdpy11.GSLrng(12345)
N = 100
L = 1000.0
ts = msprime.simulate(2*N, Ne = N, length=L, recombination_rate = 1e-3, random_seed = 42)
pop = fwdpy11.DiploidPopulation.create_from_tskit(ts)
fwdpy11.infinite_sites(rng, pop, 1e-2)
```

```{code-cell} python
ts = pop.dump_tables_to_tskit()
```

```{code-cell}
print(ts.tables.mutations[0])
```

:::{todo}
Fill in with more exposition.
:::
