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

(gvalues_fitness_vignette)=

# Genetic values - mutations with direct effects on fitness

:::{note}
The objects described here are passed to the `gvalue` parameter when initializing instances of {class}`fwdpy11.ModelParams`.
:::

In a typical population-genetic model, mutations have direct effects on fitness.
Often, this effect is referred to as `s`, or the "selection coefficient".

Once we've decided on our distributions of effect sizes, we need a way to obtain a diploid's fitness.
For these "standard" population genetic models, we will use {class}`fwdpy11.Multiplicative`.
Instances of this class tell the simulation to calculate the genetic value of an individual using a multiplicative model where the value contributed by each position with a mutation is:

| Genotype | `AA`      | `Aa`         | `aa`                        |
| -------- | --------- | ------------ | --------------------------- |
| Fitness  | {math}`1` | {math}`1+hs` | {math}`1 + scaling\times s` |

In this table:

* `A` refers to the ancestral/non-mutant allelic state
* `a` is the mutant allelic state
* `h` is the heterozygous effect of the mutant, the so-called dominance coefficient.
* `s` is the selection coefficient.
* `scaling` lets you decide between Fisher, Wright, Haldane, Kimura, etc.,
  when determining the fitness of the mutant homozygote.

The most common values for `scaling` are `1.0` or `2.0`:

```{code-cell} python
import fwdpy11

gvalue = fwdpy11.Multiplicative(scaling=1.0)
gvalue.scaling
gvalue = fwdpy11.Multiplicative(scaling=2.0)
gvalue.scaling
```

:::{note}
The `scaling` parameter interacts with the `h` parameter for a distribution of effect sizes! (See {ref}`des_vignette`.)
For example, if `scaling = 1.0`, then `h = 1.0` results in dominant mutations.
However, if `scaling = 2.0`, then `h = 1.0` gives co-dominant mutations.
In both cases, `h = 0.0` generates fully-recessive mutations.
:::
