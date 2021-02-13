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

(geneticmaps_vignette)=

# Genetic maps

:::{note}
The lists of objects described here are passed to the `recregions` parameter when initializing instances of {class}`fwdpy11.ModelParams`.
:::

:::{versionadded} 0.3.0
The methods described here replace a soon-to-be deprecated approach using {class}`fwdpy11.Region`.
:::

:::{versionchanged} 0.12.0
Added ability to restrict crossover positions to discrete values.
:::

Genetic maps are lists of instances of classes derived from the `ABC`
{class}`fwdpy11.GeneticMapUnit`.
Here `Unit` refers to an *element* of a genetic map rather than the actual units (`cM`, etc.).
Instances of these classes contain their own rates and we can mix and match regions where recombination breakpoints are Poisson and binomially distributed.

The following example uses {class}`fwdpy11.PoissonInterval` to model two genomic regions of equal length.
The first region, spanning $[0, 5)$, has a mean of $2 \times 10^{-3}$ crossovers per generation. 
The second region, spanning [5, 10), has half as many crossovers per generation.

```{code-cell} python
import fwdpy11

recregions = [
    fwdpy11.PoissonInterval(beg=0, end=5, mean=2e-3),
    fwdpy11.PoissonInterval(beg=5, end=10, mean=1e-3),
]
```

The number of breakpoints in each {math}`[beg, end)` interval is Poisson distributed with the given mean.
The position of each breakpoint is uniform (and continuous) on {math}`[beg, end)`.

These classes also allow us to specify breakpoints at a specific position with a specific probability.
The next example sets up 4 genomic regions, each 10 "units" long.  Within each region, the mean number of breakpoints (per diploid, per generation) is {math}`1e-3`.
Between each region, a single recombination occurs with probability of
one-half, meaning that each region is assorting independently (50 `cM` between each region).

```{code-cell} python
NLOCI = 4
LOCUS_LENGTH = 10
RECRATE_PER_LOCUS = 1e-3
LOCUS_BOUNDARIES = [
    (i, i + LOCUS_LENGTH) for i in range(0, NLOCI * LOCUS_LENGTH, LOCUS_LENGTH)
]
recregions = [fwdpy11.PoissonInterval(*i, RECRATE_PER_LOCUS) for i in LOCUS_BOUNDARIES]
for i in LOCUS_BOUNDARIES[:-1]:
    recregions.append(fwdpy11.BinomialPoint(i[1], 0.5))
for i in recregions:
    print(i)
```

As an aside, this example is not creating objects in order by their positions.  Such ordering is not required.

Beginning in version `0.12.0`, it is possible to restrict crossover positions to integer values.
For the examples given above, crossover positions are floating-point values sampled uniformly from {math}`[beg, end)`.
To restrict positions to integer values, we pass `discrete=True` when creating object instances:

```{code-cell} python
recregions = [
    fwdpy11.PoissonInterval(beg=0, end=5, mean=2e-3, discrete=True),
    fwdpy11.PoissonInterval(beg=5, end=10, mean=1e-3, discrete=True),
]
```

Now, breakpoints from the first region will only take on values of `0`, `1`, `2`, `3`, or `4`.

Setting `discrete=True` requires the following:

* Values for `beg` and `end` must be {class}`int`.  Thus, `1` is valid but `1.0` will raise a `TypeError`.
* `end - beg` must be `> 1`.  This requirement prevents you from using `beg=0` and `end=1`, for example, which would result in the only possible crossover position being `0`.
* You must be more careful when using `msprime` to start/finish a simulation.
  See {ref}`here <precapitation>` and {ref}`here <recapitation>` for details.

The following classes are available:

* {class}`fwdpy11.PoissonInterval`
* {class}`fwdpy11.PoissonPoint`
* {class}`fwdpy11.BinomialInterval`
* {class}`fwdpy11.BinomialPoint`
* {class}`fwdpy11.FixedCrossovers`

## General comments

* Different {math}`[beg, end)` intervals may overlap.
  The interpretation of such a setup is up to you.
* When using classes like {class}`fwdpy11.PoissonInterval`, the recombination rate that you use to construct a {class}`fwdpy11.ModelParams` instance is ignored, as the rates are stored in the individual objects.
* You do not need to specify regions with zero recombination.
  Their existence is implied given the total length of the genome being simulated ({attr}`fwdpy11.TableCollection.genome_length`).
