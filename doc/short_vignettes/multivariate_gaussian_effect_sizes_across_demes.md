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

(multivariate_gaussian_effect_sizes_across_demes_vignette)=

# The multivariate Gaussian distribution

```{code-cell}python
---
tags: ['hide-input']
---
import demes
import demesdraw
import fwdpy11
import numpy as np
```


We may model Gaussian effect sizes using the existing {class}`fwdpy11.MultivariateGaussianEffects`
in conjunction with {class}`fwdpy11.mvDES`.

At this time, it is probably best to look at an example. The following code models Gaussian stabilizing
selection on a quantitative trait.  The effects sizes within each deme are themselves given by Gaussian
distributions and there is no correlation in the effect size in the two demes.

We will simulate a demographic model of migration happening into the infinite past of two equal-sized demes:

```{code-cell}python
---
tags: ['hide-input']
---
yaml = """
description: Island model forever
time_units: generations
demes:
  - name: ancestor
    epochs:
     - start_size: 100
       end_time: 100
  - name: A
    ancestors: [ancestor]
    epochs:
     - start_size: 100
  - name: B
    ancestors: [ancestor]
    epochs:
     - start_size: 100
migrations:
  - demes: [A, B]
    rate: 0.10
"""
g = demes.loads(yaml)
model = fwdpy11.discrete_demography.from_demes(g, burnin=1)
demesdraw.tubes(g);
```

Let's set up a dictionary to hold the parameters of our model:

```{code-cell} python
pdict = {
    "nregions": [],
    "recregions": [],
    "sregions": [
        fwdpy11.mvDES(
            fwdpy11.MultivariateGaussianEffects(0, 1, 1, np.identity(3)), np.zeros(3)
        )
    ],
    "rates": (0, 2.5e-3, None),
    "demography": model,
    "simlen": model.metadata["total_simulation_length"],
    "gvalue": fwdpy11.Additive(
        ndemes=3, scaling=2,
        gvalue_to_fitness=fwdpy11.GaussianStabilizingSelection.single_trait([fwdpy11.Optimum(optimum=0.0, VS=10.0, when=0)])
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
# TODO: update this once we have a function to pull the sizes
# automatically from demes-derived models:
initial_sizes = [v for v in model.metadata["initial_sizes"].values()]
pop = fwdpy11.DiploidPopulation(initial_sizes, 1.0)
rng = fwdpy11.GSLrng(1010)
fwdpy11.evolvets(rng, pop, params, 10)
```

Let's extract the effect sizes from each deme:

```{code-cell} python
assert len(pop.tables.mutations) > 0
for i in pop.tables.mutations:
    print(pop.mutations[i.key].esizes)
```

Let's look at another example where effect sizes covary negatively across demes and raise the mutation rate a bit:

```{code-cell} python
# Effect sizes across demes will
# have a correlation coefficient of r=1/2
cor_matrix = np.array([-0.5]*9).reshape(3,3)
np.fill_diagonal(cor_matrix, np.array([1.0]*3))

# Get our covariance matrix
sd = np.array([0.1]*3)
D = np.identity(3)
np.fill_diagonal(D, sd)
vcv_matrix = np.matmul(np.matmul(D, cor_matrix), D)
pdict["sregions"] = [
    fwdpy11.mvDES(fwdpy11.MultivariateGaussianEffects(0, 1, 1, vcv_matrix), np.zeros(3))
]
params = fwdpy11.ModelParams(**pdict)
# TODO: update this once we have a function to pull the sizes
# automatically from demes-derived models:
initial_sizes = [v for v in model.metadata["initial_sizes"].values()]
pop = fwdpy11.DiploidPopulation(initial_sizes, 1.0)
fwdpy11.evolvets(rng, pop, params, 10)
for i in pop.tables.mutations:
    print(pop.mutations[i.key].esizes)
```

Now we see that the effect sizes often differ in sign between the two demes.
