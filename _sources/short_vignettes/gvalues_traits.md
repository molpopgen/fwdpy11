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

(gvalues_traits_vignette)=

# Genetic values - simulating phenotypes/traits.

:::{note}
The objects described here are passed to the `gvalue` parameter when initializing instances of {class}`fwdpy11.ModelParams`.
:::

The previous section discussed setting up a model where a mutation's effect size ({attr}`fwdpy11.Mutation.s`) directly affects individual fitness.
An alternative model is one where mutations affect some abstract "trait" or "phenotype" and a separate function maps trait values to fitness.

Let's consider the standard model of evolutionary quantitative genetics:

* Mutations have **additive** effects on trait values
* The fitness of a trait value is a quadratic function of its distance
  from an "optimum" trait value.

In `fwdpy11`, a non-mutant individual has a phenotype of `0.0`.
Trait values are additive over the values contributed by individual genotypes
according to the following table:

| Genotype    | `AA`      | `Aa`       | `aa`                    |
| ----------- | --------- | ---------- | ----------------------- |
| Trait value | {math}`0` | {math}`hs` | {math}`scaling\times s` |

(If we model multiplicative effects on a trait, a non-mutant individual still has a value of `0.0`.
The internal machinery handles this so that you don't have to worry about it.)

To specify an additive effects model of a trait under Gaussian stabilizing selection with an optimum trait value of `0.0` and (inverse) strength of stabilizing selection `VS = 1.0`, we write:

```{code-cell} python
import fwdpy11

gvalue = fwdpy11.Additive(
    scaling=2.0, gvalue_to_fitness=fwdpy11.GSS(optimum=0.0, VS=1.0)
)
```

Here, we are using a second parameter to initialize a "genetic value to fitness" map stored in an instance of {class}`fwdpy11.Additive`.
({class}`fwdpy11.Multiplicative` also supports such maps.)
See {class}`fwdpy11.GSS` for details.

We can also add Gaussian noise to an individual's trait value:

```{code-cell} python
import numpy as np

gvalue = fwdpy11.Additive(
    scaling=2.0,
    gvalue_to_fitness=fwdpy11.GSS(optimum=0.0, VS=2.0 / 3.0),
    noise=fwdpy11.GaussianNoise(mean=0.0, sd=np.sqrt(1.0 / 3.0)),
    )
```

The last example requires some explanation:

* We want `VS = 1.0`.  We can decompose `VS = VW + VE`, where `VW` and `VE` are the additive contributions of genetic and environmental effects.
* Here, the environmental effect is a Gaussian with mean zero and variance `1/3`.
  The class is parameterized with the standard deviation, however, so we need to pass on the square root.
* We then set `VS = 1 - 1/3 = 2/3` when initializing {class}`fwdpy11.GSS`.

Yes, this is a nomenclature issue!
The `VS` argument to {class}`fwdpy11.GSS` really should be called `VW` and we'll fix that in a future version and hopefully not break people's code.

In general, there's a good bit of subtlety to properly modeling quantitative traits.
The machinery described here was used in {cite}`Thornton2019-nu`. {cite}`Burger2000-ul` is an excellent technical reference on the topic.
{cite}`Walsh2018-ux` also thoroughly covers a lot of relevant material.

:::{note}

Under the hood, the `GSS` and `GSSmo` classes aren't that different.
Their multivariate analogs are rather similar, too.
Thus, we envision a future with one single `fwdpy11.GaussianStabilizingSelection` class to handle all cases.
The types discussed here would remain as simple Python wrappers so that we don't break existing simulations.

:::

For an example of another approach to modeling phenotypes often attributed to {cite}`Eyre-Walker2010-rs`, see {ref}`here <eyre-walker>`.

:::{todo}
Write (and refer to) an advanced section on pleiotropic models.
:::

## Changing the optimum phenotype during a simulation

### A sudden optimum shift

The previous example set up a model where the optimum is stable for the entire simulation.
We can parameterize a shifting optimum using {class}`fwdpy11.GSSmo`.
For example, to shift the optimum from `0.0` to `1.0` at generation `100`:

```{code-cell} python
moving_optimum = fwdpy11.GSSmo(
    [
        fwdpy11.Optimum(when=0, optimum=0.0, VS=1.0),
        fwdpy11.Optimum(when=100, optimum=1.0, VS=1.0),
    ]
)

gvalue = fwdpy11.Additive(scaling=2.0, gvalue_to_fitness=moving_optimum)
```

### Randomly-moving optimum

Since we are working in `Python`, we can take advantage of existing libraries to implement interesting models.
Let's consider the following model of a randomly moving optimum:

* There is a 1% chance each generation that the optimum shifts.
* When a shift happens, a normal deviate with mean `0.0` and variance `0.1` is added to the current optimum.
* The simulation will end at generation `1,000`.

Let's code it up:

```{code-cell} python
optima = [fwdpy11.Optimum(when=0, optimum=0.0, VS=10.0)]

last_time = 0
last_optimum = 0.0

np.random.seed(666)

while last_time < 1000:
    last_time += int(np.random.geometric(0.01, 1)[0])
    last_optimum += np.random.normal(loc=0.0, scale=np.sqrt(0.1), size=1)[0]
    if last_time < 1000:
        optima.append(fwdpy11.Optimum(when=last_time, optimum=last_optimum, VS=10.0))

random_moving_optimum = fwdpy11.GSSmo(optima)
random_moving_optimum
```

:::{note}

Note the cast to `int` when updating the time.
{class}`fwdpy11.Optimum` is very picky about its input.
It requires `int` for `when` and will raise an exception if the {attr}`numpy.int64` from {func}`numpy.random.geometric` gets passed in.

:::
