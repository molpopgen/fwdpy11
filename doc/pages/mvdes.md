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
(mvdes)=

# Different effect sizes of mutations in different demes

:::{versionadded} 0.7.0

:::

This section describes how to allow mutations to have different effect sizes in different demes.

To allow for correlated effect sizes across demes, we used {class}`fwdpy11.mvDES` to model multivariate
distributions of effect sizes.  This class inherits from {class}`fwdpy11.Sregion` (see {ref}`des_vignette`).
In general, there is no standard *general* method for drawing random deviates from arbitrary multivariate
distributions.  The approach that we take is to use a multivariate Gaussian distribution as the underlying kernel.

First, we will describe the cases where the multivariate Gaussian kernel leads to output effect size distributions
that have straightforward interpretations.  Then, we will move onto the more general case allowing you to construct
multivariate distributions from current {class}`fwdpy11.Sregion` types.

## The multivariate Gaussian

We may model Gaussian effect sizes using the existing {class}`fwdpy11.MultivariateGaussianEffects`
in conjunction with {class}`fwdpy11.mvDES`.  Using {class}`fwdpy11.MultivariateGaussianEffects` on its
own is used to model pleiotropic effects on a trait.  Here, we are re-using this type to model correlated
effect sizes across demes.

At this time, it is probably best to look at an example. The following code models Gaussian stabilizing
selection on a quantitative trait.  The effects sizes within each deme are themselves given by Gaussian
distributions and there is no correlation in the effect size in the two demes.

```{code-cell} python
import demes
import fwdpy11
import numpy as np

yaml = """
description: two demes with symmetric migration
time_units: generations
demes:
 - name: deme0
   epochs:
    - start_size: 100
 - name: deme1
   epochs:
     - start_size: 100
migrations:
 - demes: [deme0, deme1]
   rate: 0.1
"""

graph = demes.loads(yaml)
demography = fwdpy11.discrete_demography.from_demes(graph, 1)

pdict = {
    "nregions": [],
    "recregions": [],
    "sregions": [
        fwdpy11.mvDES(
            fwdpy11.MultivariateGaussianEffects(0, 1, 1, np.identity(2)), np.zeros(2)
        )
    ],
    "rates": (0, 1e-3, None),
    "demography": demography,
    "simlen": 100,
    "gvalue": fwdpy11.Additive(
        ndemes=2, scaling=2, gvalue_to_fitness=fwdpy11.GSS(optimum=0.0, VS=10.0)
    ),
}
```

Most of the above is standard.  Let's dissect the new bits:

* An instance of {class}`fwdpy11.mvDES` is our only region with selected mutations.
* This instance holds an instance of {class}`fwdpy11.MultivariateGaussianEffects`
  that puts mutations on the interval {math}`[0, 1)` with weight 1 and an identity
  matrix specifies the correlation in effect sizes between demes 0 and 1.  The
  identity matrix has the value zero for all off-diagonal elements, meaning
  no covariance in effect sizes across demes.
* The final constructor argument specifies the mean of each marginal Gaussian
  distribution. The means are both zero.
* Our genetic value type accepts an `ndemes` parameter, telling it that it has
  to look for deme-specific effect sizes.  This value must be set to the maximum
  number of demes that will exist during a simulation.

Let's evolve the model now:

```{code-cell} python
params = fwdpy11.ModelParams(**pdict)
pop = fwdpy11.DiploidPopulation([100, 100], 1.0)
rng = fwdpy11.GSLrng(1010)
fwdpy11.evolvets(rng, pop, params, 10)
```

Let's extract the effect sizes from each deme:

```{code-cell} python
for i in pop.tables.mutations:
    print(pop.mutations[i.key].esizes)
```

Let's look at another example where effect sizes covary negatively across demes and raise the mutation rate a bit:

```{code-cell} python
vcv = np.array([1.0, -0.99, -0.99, 1.0]).reshape((2, 2))
pdict["sregions"] = [
    fwdpy11.mvDES(fwdpy11.MultivariateGaussianEffects(0, 1, 1, vcv), np.zeros(2))
]
pdict["rates"] = (0, 5e-3, None)
params = fwdpy11.ModelParams(**pdict)
pop = fwdpy11.DiploidPopulation([100, 100], 1.0)
fwdpy11.evolvets(rng, pop, params, 10)
for i in pop.tables.mutations:
    print(pop.mutations[i.key].esizes)
```

Now we see that the effect sizes often differ in sign between the two demes.

## The multivariate lognormal

If {math}`X` is a multivariate Gaussian distribution, {math}`N(\mathbf{\mu}, \mathbf{\sum})`, where {math}`\mathbf{\mu}` is a vector of mean values and
{math}`\mathbf{\sum}` is the covariance matrix, then {math}`Y = e^X` is a
multivariate lognormal random variable with mean {math}`E[Y]_i = e^{\mu_i + \frac{1}{2}\sum_{ii}}` and covariance matrix {math}`Var[Y]_{i,j} = e^{\mu_i + \mu_j + \frac{1}{2}(\sum_{ii} + \sum_{jj})}(e^{\sum_{ij}}-1)`.

To specify a multivariate lognormal distribution of effect sizes, we use
the static class method {func}`fwdpy11.LogNormalS.mv`.  The following code
constructs a distribution of effect sizes such that `-2Ns` (where `N` is the
size of a single deme) is a multivariate lognormal with means zero and an
identity matrix as a covariance matrix used to specify the multivate
Gaussian kernel.

```{code-cell} python
mvdes = fwdpy11.mvDES(
    fwdpy11.LogNormalS.mv(0, 1, 1, scaling=-200), np.zeros(2), np.identity(2)
)
```

:::{note}

The lognormal distribution returns deviates {math}`> 0`.
To model deleterious mutations/effect sizes < 0, use the
`scaling` parameter with a negative value like we just did!

:::

Let's put it in a simulation and run it:

```{code-cell} python
pdict = {
    "nregions": [],
    "recregions": [],
    "sregions": [mvdes],
    "rates": (0, 1e-3, None),
    "demography": demography,
    "simlen": 100,
    "gvalue": fwdpy11.Multiplicative(ndemes=2, scaling=2),
}
params = fwdpy11.ModelParams(**pdict)
pop = fwdpy11.DiploidPopulation([100, 100], 1.0)
fwdpy11.evolvets(rng, pop, params, 10)
for i in pop.tables.mutations:
    print(pop.mutations[i.key].esizes)
```

## "Custom" multivariate distributions

The previous two sections cover cases where the methods for generating
deviates from a multivariate distribution are straightforward and agreed
upon.

In order to simulate multivariate distributions of effect sizes based on
{class}`fwdpy11.Sregion` types, we follow a fairly intuitive approach
described in {cite}`Xue-Kun_Song2000-qn`.  Briefly, the multivariate Gaussian kernel is
used to produce deviates.  Then, the quantiles from the cummulative distribution
of each marginal Gaussian are used to generate a deviate from the desired output distribution of interest.

For a simulation with `n` populations we need:

* A {class}`list` of `n` {class}`fwdpy11.Sregion` objects
* An array of `n` means for the multivariate Gaussian
* An `n-by-n` covariance matrix for the multivariate
  Gaussian

The following generates exponentially distributed effect sizes in each deme
with a high correlation across demes:

```{code-cell} python
mvdes = fwdpy11.mvDES(
    [fwdpy11.ExpS(0, 1, 1, -0.5)] * 2,
    np.zeros(2),
    np.matrix([1, 0.9, 0.9, 1]).reshape((2, 2)),
)
pdict["sregions"] = [mvdes]
params = fwdpy11.ModelParams(**pdict)
pop = fwdpy11.DiploidPopulation([100, 100], 1.0)
fwdpy11.evolvets(rng, pop, params, 10)
for i in pop.tables.mutations:
    print(pop.mutations[i.key].esizes)
```

We can mix and match our distributions.  Here, the distribution of effect
sizes in deme 0 is exponential and the distribution in deme 1 is gamma.  The
two distributions have means with opposite signs and the magnitudes of the
marginal deviates negatively covary:

```{code-cell} python
mvdes = fwdpy11.mvDES(
    [fwdpy11.ExpS(0, 1, 1, -0.5), fwdpy11.GammaS(0, 1, 1, mean=0.1, shape_parameter=1)],
    np.zeros(2),
    np.matrix([1, -0.9, -0.9, 1]).reshape((2, 2)),
)
pdict["sregions"] = [mvdes]
params = fwdpy11.ModelParams(**pdict)
pop = fwdpy11.DiploidPopulation([100, 100], 1.0)
fwdpy11.evolvets(rng, pop, params, 10)
for i in pop.tables.mutations:
    print(pop.mutations[i.key].esizes)
```

The type {class}`fwdpy11.ConstantS` has intuitive behavior:

```{code-cell} python
mvdes = fwdpy11.mvDES(
    [fwdpy11.ExpS(0, 1, 1, -0.5), fwdpy11.ConstantS(0, 1, 1, -0.1)],
    np.zeros(2),
    np.matrix([1, -0.9, -0.9, 1]).reshape((2, 2)),
)
pdict["rates"] = (0, 5e-3, None)
pdict["sregions"] = [mvdes]
params = fwdpy11.ModelParams(**pdict)
pop = fwdpy11.DiploidPopulation([100, 100], 1.0)
rng = fwdpy11.GSLrng(1010)
fwdpy11.evolvets(rng, pop, params, 10)
for i in pop.tables.mutations:
    print(pop.mutations[i.key].esizes)
```

## Recipes

### Different signs in different demes

Consider two demes.  You want any beneficial mutation in one deme to
be deleterious in the other and vice-versa.

For the multivariate Gaussian, use the covariance matrix as done above.  Note
that this approach only generates a *tendency* to different signs in different demes.

With the multivariate lognormal, the best we can do is to use negative
correlations such that deleterious mutations in deme 0 are less deleterious in deme 1, etc.:

```{code-cell} python
sregions = [
    fwdpy11.mvDES(
        fwdpy11.LogNormalS.mv(0, 1, 1, scaling=-200),
        np.zeros(2),
        np.matrix([1, -0.99, -0.99, 1]).reshape((2, 2)),
    )
]
sregions.append(
    fwdpy11.mvDES(
        fwdpy11.LogNormalS.mv(0, 1, 1, scaling=200),
        np.zeros(2),
        np.matrix([1, -0.99, -0.99, 1]).reshape((2, 2)),
    )
)
pdict["sregions"] = sregions
params = fwdpy11.ModelParams(**pdict)
pop = fwdpy11.DiploidPopulation([100, 100], 1.0)
rng = fwdpy11.GSLrng(1010)
fwdpy11.evolvets(rng, pop, params, 10)
for i in pop.tables.mutations:
    print(pop.mutations[i.key].esizes)
```

In the output, we see that an effect size in deme `i` has a corresponding effect size in deme `j` that is a about an order of magnitude smaller in absolute value.

For the general approach, simply create a {class}`list` of objects with the desired mean (or constant) effect sizes.  For example:

```{code-cell} python
sregions = [
    fwdpy11.mvDES(
        [fwdpy11.ExpS(0, 1, 1, -0.5), fwdpy11.ExpS(0, 1, 1, 0.1)],
        np.zeros(2),
        np.identity(2),
    )
]
sregions.append(
    fwdpy11.mvDES(
        [fwdpy11.ExpS(0, 1, 1, 0.1), fwdpy11.ExpS(0, 1, 1, -0.5)],
        np.zeros(2),
        np.identity(2),
    )
)
pdict["sregions"] = sregions
params = fwdpy11.ModelParams(**pdict)
pop = fwdpy11.DiploidPopulation([100, 100], 1.0)
rng = fwdpy11.GSLrng(1010)
fwdpy11.evolvets(rng, pop, params, 10)
for i in pop.tables.mutations:
    print(pop.mutations[i.key].esizes)
```

### Polygenic traits, multiple demes, correlated effect sizes, and different optima

See {ref}`GSSDivergentOptima`.


