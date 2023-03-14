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

(workingexample_fitness_vignette)=

# Working example - parameters for simulating direct effects on fitness

```{code-cell}
import fwdpy11

N = 1000
mu_neutral = 0.0
mu_selected = 1e-3
rho = 1000.
alpha = -50

p = {
    "nregions": [],
    "gvalue": fwdpy11.Multiplicative(2.0),
    "sregions": [fwdpy11.ExpS(0, 1., 1, mean=alpha/(2 * N))],
    "recregions": [fwdpy11.PoissonInterval(0, 1., rho / float(4 * N))],
    "rates": (mu_neutral, mu_selected, None),
    "prune_selected": True,
    "demography": None,
    "simlen": 10 * N,
}
params = fwdpy11.ModelParams(**p)
```

Let's take a look at what he have:

```{code-cell}
print(params.asblack())
```

Let's explain a few new things:

* The `rate` {class}`tuple` contains the neutral mutation rate, selected mutation rate, and recombination rate, respectively.
  The recombination rate is `None` because the rates are all contained in `recregions`.
  (See {ref}`geneticmaps_vignette`.)
  This parameterization of `rates` will soon be deprecated and the need to specify a recombination rate of `None` will go away.
* The `demography` field is something we've not encountered yet.
  If nothing is provided, an empty instance of {class}`fwdpy11.ForwardDemesGraph` gets initialized.
  This means "no demographic events".
  Demographic events are the subject of {ref}`this page <demes_vignette>`.
* `simlen` is how many generations to evolve the population.
* `prune_selected` is a boolean.
  If `True`, then mutations that are fixed in the alive individuals, and not present at all in ancient samples, will be removed from the simulation.
  This removal is an efficiency thing--we don't need them for multiplicative models.
  In the future, we hope for a better method for handling fixations that allows for removal, even when present in ancient samples, during the simulation, then add them back in to the tables before returning from the simulation.
  
:::{warning}
Do not set `prune_selected = True` for models of additive effects.
Doing so breaks our requirement that the rank order of relative fitness is preserved up to a multiplicative constant.
:::

Our next vignette shows us some tricks that can be useful for modifying our `params` objects.
