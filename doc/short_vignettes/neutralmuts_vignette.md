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

(neutralmuts_vignette)=

# Neutral mutations during a simulation

:::{note}
The lists of objects described here are passed to the `nregions` and `rates` parameters when initializing instances of {class}`fwdpy11.ModelParams`.
:::

Normally, there is no need to simulate neutral mutations during a simulation.
Because we record the ancestry using tree sequence recording, it is much more efficient in general to add neutral mutations when the simulation is complete.
See {ref}`here <neutralmutsafter_vignette>`.

The reasons to simulate neutral variants during a simulation include:

* Wanting to track their frequencies over time.
* Wanting neutral mutation rates to vary along the genome.
  (This is only necessary due to a current limitation of {func}`fwdpy11.infinite_sites`.)

The procedure is relatively simple.
Variation in the neutral mutation rate is specified by different weights in instances of {class}`fwdpy11.Region`.

For example, to have two regions of the same size, but the second having twice as many mutations as the first:

```{code-cell}
import fwdpy11

nregions = [fwdpy11.Region(beg=0.0, end=5.0, weight=1.0),
            fwdpy11.Region(beg=5.0, end=10.0, weight=2.0)
           ]
```

To specify the total mutation rate to neutral mutations, pass a non-negative {class}`float` as the first element of `rates`.
(See {class}`fwdpy11.ModelParams`.)

To see neutral mutations during a simulation in action, check out {ref}`this vignette <freqtrackpy_neutral>`.
