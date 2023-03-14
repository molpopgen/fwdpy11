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

(neutralmutsafter_vignette)=

# Adding neutral mutations to a population with existing ancestry

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

We can add mutations to the tables via {func}`fwdpy11.infinite_sites`:

```{code-cell}
nmuts = fwdpy11.infinite_sites(rng, pop, 1e-3)
print(f"{nmuts} mutations added")
```

We can check that the tables do indeed have these variants:

```{code-cell}
neut = 0
non_neut = 0
for m in pop.tables.mutations:
    if m.neutral is True:
        neut += 1
    else:
        non_neut += 1
print(f"There are {neut} neutral and {non_neut} non-neutral mutations")
```

The number of genomes containing each mutation get recorded:

```{code-cell}
for m in pop.tables.mutations:
    if m.neutral is True:
        print(pop.mcounts[m.key])
```

## Limitations

* Only the infinitely-many sites model is implemented here.
  For more general models, you can transfer the data to `tskit` format and use that library to add mutations.
  See {ref}`tskitconvert_vignette`.
* Currently, we only support a constant mutation rate across the genome.
  We will generalize this in a future release.
  In the mean time, you can add neutral mutations during the simulation as described {ref}`here <neutralmuts_vignette>`.

