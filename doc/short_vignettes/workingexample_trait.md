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

(workingexample_trait_vignette)=

# Working example - parameters for the simulation of a trait

```{code-cell}
import fwdpy11

N = 1000
rho = 1000.0
mu_neutral = 0.0
mu_selected = 1e-3

GSSmo = fwdpy11.GSSmo(
    [
        fwdpy11.Optimum(when=0, optimum=0.0, VS=1.0),
        fwdpy11.Optimum(when=10 * N, optimum=1.0, VS=1.0),
    ]
)

p = {
    "nregions": [],
    "gvalue": fwdpy11.Additive(2.0, GSSmo),
    "sregions": [fwdpy11.GaussianS(0, 1., 1, 0.1)],
    "recregions": [fwdpy11.PoissonInterval(0, 1., rho / float(4 * N))],
    "rates": (mu_neutral, mu_selected, None),
    # Keep mutations at frequency 1 in the pop if they affect fitness.
    "prune_selected": False,
    "demography": fwdpy11.DiscreteDemography(),
    "simlen": 10 * N + 200,
}
params = fwdpy11.ModelParams(**p)
```

```{code-cell}
print(params.asblack())
```

## Modifying model parameters

Instances of {class}`fwdpy11.ModelParams` are immutable.
To change them, either:

* Modify the original `dict`
* Or, do a round trip through a temporary `dict`.

The first option is straightforward.
Let's see the second in action.
We'll change the simulation length.

```{code-cell}
temp_dict = params.asdict()
temp_dict["simlen"] *= 2
params = fwdpy11.ModelParams(**temp_dict)
print(params.simlen)
