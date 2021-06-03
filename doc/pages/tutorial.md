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

(tutorial)=

# Tutorial

The tutorial is organized around a series of tasks that one
needs to complete in order to run a simulation.

Each section covers a topic relevant to parameterizing a simulation.
The variable names that store the types are the same as `kwargs`
used to initialize the final model parameters class.

Along the way, you will see links to more detailed and/or more advanced
sections of the documentation.

For a first-time reader, it is probably best to go top-to-bottom through
this tutorial.

## Initializing a population

A simulation modifies the contents of instances of {class}`fwdpy11.DiploidPopulation`.
We have two basic ways to initialize a population:

* with a list of deme sizes and a genome length
* from a tree sequence in `tskit` format

This section deals with the former.  See {ref}`here <pop-from-tskit>`
for how to initialize from {class}`tskit.TreeSequence` objects.

To initialize a population of `N` diploids from a single deme and a given
genome length:

```{code-cell} python
import fwdpy11
import numpy as np

pop = fwdpy11.DiploidPopulation(N=1000, length=100.0)
pop.N
pop.tables.genome_length
```

The genome length is stored as an attribute of the *table collection*
that is used to record the ancestry of the simulated population.
The tables themselves are instances of {class}`fwdpy11.TableCollection`.
See {ref}`here <tsoverview>` for a brief introduction to tables as well
as {cite}`Kelleher2016-cb` and {cite}`Kelleher2018-mc`.  See {ref}`here <ts-data-types>`
for a more detailed discussion of relevant data structures.

To initialize a population with individuals in multiple demes:

```{code-cell} python
pop = fwdpy11.DiploidPopulation(N=[500, 500], length=100.0)
pop.N
pop.tables.genome_length
```

### Looking at individual metadata

Our population size is unchanged, but our population now has 500 individuals
in each of two demes.  We can verify this by looking at the `metadata` associated
with each individual.  Individual metadata are instances of
{class}`fwdpy11.DiploidMetadata` but they are also registered as a {class}`numpy.dtype`,
meaning that we can look at the metadata as a {class}`numpy.recarray`:

```{code-cell} python
md = np.array(pop.diploid_metadata, copy=False)
```

The field names of our array are the same as the attribute names of
{class}`fwdpy11.DiploidMetadata`:

```{code-cell} python
md.dtype
```

We can easily confirm the number of individuals in each deme using
{func}`numpy.unique`:

```{code-cell} python
np.unique(md["deme"], return_counts=True)
```

We see that the deme labels are `0` and `1` and that each label
was found 500 times.  (The first 500 individuals are in deme `0`,
followed by 500 in deme `1`.)

See {ref}`here <tskit_metadata_vignette>` for more on individual metadata.

(mutationregions)=

## Defining distributions of mutation effect sizes

One of the main reasons to perform forward simulations is to be able
to model mutations affecting individual fitness. To do so, we need
to specify both mutation rates and the resulting effect sizes.

`fwdpy11` works by specifying an overall mutation rate to variants
affecting fitness (see {ref}`here <model-params>`).  Given that
a mutation occurs, we need to specify its "effect size".

`fwdpy11` chooses the effect size of a new mutation by first
determining what *region* is mutated and then generating a mutation
from the distribution of effect size associated with that region.

Each region is represented by instances of classes derived from the
ABC {class}`fwdpy11.Sregion`.  Each instance is associated with a *weight*.
These weights establish the relative probability that a mutation comes
from a given region.  Thus, given an overall mutation
rate to non-neutral variants, instances of "`sregions`" are used to set up
a multinomial distribution for generating new mutations.

The following sets up a model where mutations have a constant effect size
({math}`s=-0.01`), dominance {math}`h=0.25`, and occur uniformly on
the interval {math}`[0, 1)`:

```{code-cell} python
sregions = [fwdpy11.ConstantS(beg=0.0, end=1.0, weight=1.0, s=-0.01, h=0.25)]
```

The previous example uses argument names for clarity, and the following is equivalent,
with the `int` values getting converted to `float` automatically:

```{code-cell} python
sregions = [fwdpy11.ConstantS(0, 1, 1, -0.01, 0.25)]
sregions[0]
```

Note that the constructor parameters for these classes often have default
values--see the specific class documentation for details.

In some scenarios, it is useful to think about the distribution of effect sizes
as scaled with respect to the population size.  For example, selection coefficients
may be exponentially-distributed with a mean of {math}`2Ns`.  To do this in
`fwdpy11`:

```{code-cell} python
# ALPHA = 2Ns
MEAN_ALPHA = -10
N = 1000
sregions = [fwdpy11.ExpS(0, 1, 1, MEAN_ALPHA, scaling=2 * N)]
sregions[0]
```

### Region weighting

When multiple "sregion" objects are used, the default behavior is to multiply
the input `weight` by `end-beg`:

```{code-cell} python
sregions = [
   fwdpy11.ExpS(beg=0.0, end=1.0, weight=1.0, mean=-0.2),
   fwdpy11.ConstantS(beg=1.0, end=3.0, weight=1.0, s=-0.1),
]
sregions
```

Here, the input `weight` is interpreted to mean the weight "per site" is constant.
In this example, twice as many mutations will have positions in {math}`[1, 3)` as from {math}`[0, 1)`.

To change the default behavior, one can prevent the coupling between input `weight` and region length:

```{code-cell} python
sregions = [
   fwdpy11.ExpS(beg=0.0, end=1.0, weight=1.0, mean=-0.2, coupled=False),
   fwdpy11.ConstantS(beg=1.0, end=3.0, weight=1.0, s=-0.1, coupled=False),
]
sregions
```

The absolute values of the `weight` parameters themselves is irrelevant.
The only thing that matters is the *relative* values from region to region.
Simulations based on the above examples would give the same results if the `weight` were 42 or 73.06.
Therefore, we can recreate our first example with code like the following:

```{code-cell} python
sregions = [
   fwdpy11.ExpS(beg=0.0, end=1.0, weight=56.0, mean=-0.2, coupled=False),
   fwdpy11.ConstantS(beg=1.0, end=3.0, weight=112.0, s=-0.1, coupled=False),
]
sregions
```

In the above example, twice as many mutations occur in the second region
because the weights have relative values of 2:1.

:::{note}

Different regions are allowed to overlap, allowing the simulation of
concepts like "coding regions" where the DFE are a weighted mixture
from multiple distributions, etc.

:::

### Setting the dominance of new mutations

The dominance of a new mutation is set by the `h` parameter during
initialization:

```{code-cell} python
fwdpy11.ExpS(beg=0.0, end=1.0, weight=1.0, mean=-0.2, h=1.0)
fwdpy11.ExpS(beg=0.0, end=1.0, weight=1.0, mean=-0.2, h=0.0)
```

### Built-in distributions of effect sizes

* {class}`fwdpy11.ConstantS`
* {class}`fwdpy11.UniformS`
* {class}`fwdpy11.ExpS`
* {class}`fwdpy11.GammaS`
* {class}`fwdpy11.GaussianS`
* {class}`fwdpy11.MultivariateGaussianEffects`
* {class}`fwdpy11.LogNormalS`
* {class}`fwdpy11.DiscreteDESD`

### Labelling mutations from different regions

It may be of use to know what region a mutation came from.  To do
so, give a nonzero value to the `label` argument:

```{code-cell} python
fwdpy11.ConstantS(beg=0.0, end=1.0, weight=1.0, s=1e-3, label=1)
```

At the end of the simulation, mutations from this region will have
the `label` value stored in the attribute {attr}`fwdpy11.Mutation.label`.

The value of `label` must fit into a 16-bit unsigned integer,
*e.g.*, {attr}`numpy.uint16`. Larger values, or negative values, will result
in exceptions.  The following example tries to use a value one larger than
the maximum allowed:

```{code-cell} python
---
tags: [raises-exception]
---
fwdpy11.ConstantS(
    beg=0.0, end=1.0, weight=1.0, s=1e-3, label=np.iinfo(np.uint16).max + 1
)
```

(geneticmaps)=

## Modeling recombination

Recombination rates are allowed to vary along genomes in a discrete fashion.  `fwdpy11`
provides two complementary methods for setting up such variation.

(recregions)=

### Method 1: regions and weights

The method described in this section works in combination with a total overall recombination
rate.  This rate is the mean of a Poisson distribution and the intervals where recombination
breakpoints happen are chosen based on their relative weights.  The regions are instances
of {class}`fwdpy11.Region`.  A region represents a continuous, half-open interval within which
crossover positions are uniformly distributed.

By way of example, say we want the following genetic map:

* The total recombination rate per diploid is {math}`1e-3`, which is the mean of a Poisson process.
* Our genome is continuous on {math}`[0,10)`.
* The recombination rate is twice as high in one part of the "genome" than in the other.

To initialize a region object, the following parameters may be used:

* `beg`, the start of the region
* `end`, the end of the region
* `weight`, the weight assigned to the region
* `coupled` is a `bool` and determines how the weights are handled internally (see below).
* `label` is an integer that defaults to `0` and is not relevant to recombination.

The first three parameters are required.  A valid region has {math}`beg \geq 0`,
{math}`end > beg` and {math}`weight >= 0` and defines a half-open interval {math}`[beg, end)`.

:::{note}

A `weight` of `0` is the same as simply not defining a region! There is
no requirement that all genetic map elements cover the entire genome. We allow
zero-weight regions for those who think that it is cleaner/more explicit to write
them down.

:::

By default, `coupled=True`, which means that the *total* weight assigned to a region
will be {math}`weight\times (end-beg)`.  It is helpful to view {math}`weight` as
the "rate per unit" and {math}`end-beg` as the number of units in the region. (For example,
"unit" could refer to base pairs, but it need not.)

There are two ways to set this model up.  The first is arguably the most intuitive, which is to make
one region twice as long as the other:

```{code-cell} python
import fwdpy11

recrate = 1e-3
recregions = [
    fwdpy11.Region(beg=0.0, end=10.0 / 3.0, weight=1.0),
    fwdpy11.Region(beg=10.0 / 3.0, end=10, weight=1.0),
]
for r in recregions:
    print(r)
```

In the output, you see that `coupled=True`, which means that the simulation's
back-end will assign twice as many crossovers to the second region as to the first.

In words, the `recregions` list and `recrate` value mean the following:

* The number of crossovers per diploid is Poisson distributed with mean
  `recrate`, or `0.001`. See {ref}`here <model-params-rate-details>` for how to
  send the `recrate` to a simulation.
* Each crossover breakpoint has a `1/3` chance of being uniformly
  distributed in {math}`[0, 10/3)` and a `2/3` chance of being
  uniformly distributed in {math}`[10/3, 10)`.

A more abstract approach relies on setting `coupled=False`, which means
that the "raw" weights that you input are the exact values used internally:

```{code-cell} python
recregions = [
    fwdpy11.Region(beg=0, end=5, weight=1, coupled=False),
    fwdpy11.Region(beg=5, end=10, weight=2, coupled=False),
]
for r in recregions:
    print(r)
```

Now, the `weight` arguments are treated as *absolute*, or exactly `1` and `2`, respectively.

In words, what we have is:

* The number of breakpoints per diploid is Poisson distributed with mean {math}`1e-3`
* For each breakpoint, its position is uniform on {math}`[0, 5)` with probability {math}`2/(2+1)`, or
  it is uniform on {math}`[5, 10)` with probability {math}`1/(2+1)`.

In essence, instances of {class}`fwdpy11.Region` parameterize a multinomial distribution that is used to
choose the ranges within which breakpoints are uniformly-distributed.  A limitation of this approach
is that we cannot model discrete jumps in genetic maps, such as those between chromosomes.

(geneticmapunit)=

### Method 2: using "genetic map" classes

:::{versionadded} 0.3.0

:::

An alternate approach uses instances of classes derived from the `ABC`
{class}`fwdpy11.GeneticMapUnit`. Here `Unit` refers to an *element* of
a genetic map rather than the actual units (`cM`, etc.).  Instances of
these classes contain their own rates and we can mix and match regions
where recombination breakpoints are Poisson and binomially distributed.

Let's revisit the example from the previous section.  This time, we will
use {class}`fwdpy11.PoissonInterval`:

```{code-cell} python
recregions = [
    fwdpy11.PoissonInterval(beg=0, end=5, mean=2e-3 / 3),
    fwdpy11.PoissonInterval(beg=5, end=10, mean=1e-3 / 3),
]
```

The number of breakpoints in each {math}`[beg, end)` interval is Poisson distributed
with the given mean. The position of each breakpoint is uniform on {math}`[beg, end)`.

These classes also allow us to specify breakpoints at a specific position with a specific probability.
The next example sets up 4 genomic regions, each 10 "units" long.  Within each region,
the mean number of breakpoints (per diploid, per generation) is {math}`1e-3`.
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
    fwdpy11.PoissonInterval(beg=0, end=5, mean=2e-3 / 3, discrete=True),
    fwdpy11.PoissonInterval(beg=5, end=10, mean=1e-3 / 3, discrete=True),
]
```

Now, breakpoints from the first region will only take on values of `0`, `1`, `2`, `3`, or `4`.

Setting `discrete=True` requires the following:

* Values for `beg` and `end` must be {class}`int`.  Thus, `1` is valid but `1.0` will raise a `TypeError`.
* `end - beg` must be `> 1`.  This requirement prevents you from using `beg=0` and `end=1`, for example,
  which would result in the only possible crossover position being `0`.
* You must be more careful when using `msprime` to start/finish a simulation.
  See {ref}`here <precapitation>` and {ref}`here <recapitation>` for details.

The following classes are available:

* {class}`fwdpy11.PoissonInterval`
* {class}`fwdpy11.PoissonPoint`
* {class}`fwdpy11.BinomialInterval`
* {class}`fwdpy11.BinomialPoint`
* {class}`fwdpy11.FixedCrossovers`

### General comments

* Different {math}`[beg, end)` intervals may overlap.  The interpretation of such a setup is your problem.
* The first method, based on {class}`fwdpy11.Region` is slightly faster, but less flexible.  More on the flexibility
  below.
* When using classes like {class}`fwdpy11.PoissonInterval`, the recombination rate that you use to construct a
  {class}`fwdpy11.ModelParams` instance is ignored, as the rates are stored in the individual objects.
* You do not need to specify regions with zero recombination. Their existence is implied given the total
  length of the genome being simulated ({attr}`fwdpy11.TableCollection.genome_length`).

:::{note}

Adding neutral mutations to the tables with {func}`fwdpy11.infinite_sites` will place
neutral variants in the non-recombining regions.

:::

(genetic-values)=

## Modeling mutations having direct effects on fitness

In a typical population-genetic model, mutations have direct effects on fitness.
Often, this effect is referred to as `s`, or the "selection coefficient".

Once we've decided on our distributions of effect sizes, we need a way to obtain
a diploid'd fitness.  For these "standard" population genetic models, we will use
{class}`fwdpy11.Multiplicative`.  Instances of this class tell the simulation
to calculate the genetic value of an individual using a multiplicative model where the
value contributed by each position with a mutation is:

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
gvalue = fwdpy11.Multiplicative(scaling=1.0)
gvalue.scaling
gvalue = fwdpy11.Multiplicative(scaling=2.0)
gvalue.scaling
```

:::{note}

The `scaling` parameter interacts with the `h` parameter
for a distribution of effect sizes! (See {ref}`mutationregions`.)
For example, if `scaling = 1.0`, then `h = 1.0` results
in dominant mutations.  However, if `scaling = 2.0`, then
`h = 1.0` gives co-dominant mutations.  In both cases,
`h = 0.0` generates fully-recessive mutations.

:::

(modeling-quant-traits)=

## Modeling mutations affecting phenotypes

The previous section discussed setting up a model where a mutation's
effect size ({attr}`fwdpy11.Mutation.s`) directly affects individual fitness.
An alternative model is one where mutations affect some abstract "trait"
or "phenotype" and a separate function maps trait values to fitness.

Let's consider the standard model of evolutionary quantitative genetics:

* Mutations have **additive** effects on trait values
* The fitness of a trait value is a quadratic function of its distance
  from an "optimum" trait value.

In `fwdpy11`, a non-mutant individual has a phenotype of `0.0`.  Trait
values are additive over the values contributed by individual genotypes
according to the following table:

| Genotype    | `AA`      | `Aa`       | `aa`                    |
| ----------- | --------- | ---------- | ----------------------- |
| Trait value | {math}`0` | {math}`hs` | {math}`scaling\times s` |

(If we model multiplicative effects on a trait, a non-mutant individual
still has a value of `0.0`. The internal machinery handles this so
that you don't have to worry about it.)

To specify an additive effects model of a trait under Gaussian
stabilizing selection with an optimum trait value of `0.0` and
(inverse) strength of stabilizing selection `VS = 1.0`, we write:

```{code-cell} python
gvalue = fwdpy11.Additive(
    scaling=2.0, gvalue_to_fitness=fwdpy11.GSS(optimum=0.0, VS=1.0)
)
```

Here, we are using a second parameter to initialize a "genetic value to fitness"
map stored in an instance of {class}`fwdpy11.Additive`.
({class}`fwdpy11.Multiplicative` also supports such maps.)
See {class}`fwdpy11.GSS` for details.

We can also add Gaussian noise to an individual's trait value:

```{code-cell} python
gvalue = fwdpy11.Additive(
    scaling=2.0,
    gvalue_to_fitness=fwdpy11.GSS(optimum=0.0, VS=2.0 / 3.0),
    noise=fwdpy11.GaussianNoise(mean=0.0, sd=np.sqrt(1.0 / 3.0)),
    )
```

The last example requires some explanation:

* We want `VS = 1.0`.  We can decompose `VS = VW + VE`, where `VW` and
  `VE` are the additive contributions of genetic and environmental effects.
* Here, the environmental effect is a Gaussian with mean zero and variance
  `1/3`.  The class is parameterized with the standard deviation, however,
  so we need to pass on the square root.
* We then set `VS = 1 - 1/3 = 2/3` when initializing {class}`fwdpy11.GSS`.

Yes, this is a nomenclature issue!  The `VS` argument to {class}`fwdpy11.GSS`
really should be called `VW` and we'll fix that in a future version and hopefully
not break people's code.

In general, there's a good bit of subtlety to properly modeling quantitative traits.
The machinery described here was used in {cite}`Thornton2019-nu`. {cite}`Burger2000-ul` is an excellent
technical reference on the topic. {cite}`Walsh2018-ux` also thoroughly covers a lot of
relevant material.

:::{note}

Under the hood, the `GSS` and `GSSmo` classes aren't that different.
Their multivariate analogs are rather similar, too.  Thus, we envision
a future with one single `fwdpy11.GaussianStabilizingSelection` class
to handle all cases.  The types discussed here would remain as simple
Python wrappers so that we don't break existing simulations.

:::

For an example of another approach to modeling phenotypes often
attributed to {cite}`Eyre-Walker2010-rs`, see {ref}`here <eyre-walker>`.

:::{todo}
    Write (and refer to) an advanced section on pleiotropic models.
:::

### Changing the optimum phenotype during a simulation

The previous example set up a model where the optimum is stable for
the entire simulation.  We can parameterize a shifting optimum
using {class}`fwdpy11.GSSmo`.  For example, to shift the optimum from `0.0`
to `1.0` at generation `100`:

```{code-cell} python
moving_optimum = fwdpy11.GSSmo(
    [
        fwdpy11.Optimum(when=0, optimum=0.0, VS=1.0),
        fwdpy11.Optimum(when=100, optimum=1.0, VS=1.0),
    ]
)

gvalue = fwdpy11.Additive(scaling=2.0, gvalue_to_fitness=moving_optimum)
```

Since we are working in `Python`, we can take advantage of existing libraries to
implement interesting models.  Let's consider the following model of a randomly
moving optimum:

* There is a 1% chance each generation that the optimum shifts.
* When a shift happens, a normal deviate with mean `0.0` and variance
  `0.1` is added to the current optimum.
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

Note the cast to `int` when updating the time.  {class}`fwdpy11.Optimum`
is very picky about its input. It requires `int` for `when` and will
raise an exception if the {attr}`numpy.int64` from {func}`numpy.random.geometric`
gets passed in.

:::

## Adding demographic events involving discrete demes

`fwdpy11` has a flexible interface for demographic models involving multiple
discrete demes.  A full overview of the `API` is given {ref}`here <softselection>`.
This section gives a cursory introduction.

Consider the following verbal description of a model:

* There is a single ancestral population of `N = 100` diploids.
* 100 generations into this population's future, it splits into two equal-sized
  demes
* The migration rate between each deme is `1e-3`.

To set this model up, first initialize the population to the ancestral state:

```{code-cell} python
pop = fwdpy11.DiploidPopulation(100, 1.0)
```

This model has migration, so we need a migration matrix. The rows of a migration
matrix are the **destination** demes and the columns are the **source** demes.
A migration matrix can be interpreted as the fraction of migrants each generation
from each source deme. This definition implies that each row must sum to `1.0`.

This is a 2-deme model, so we need a `2x2` matrix.  Initially, there is only the
single ancestral deme, and therefore 100% of its ancestry each generation is
from itself. Thus, our initial migration matrix looks like:

```{code-cell} python
migmatrix = np.zeros(4).reshape(2, 2)
migmatrix[0, 0] = 1.0
migmatrix
```

At generation `100`, we move half of the ancestral population (deme `0`) to a new
deme `1`:

```{code-cell} python
mass_migrations = [
    fwdpy11.move_individuals(when=100, source=0, destination=1, fraction=0.5)
]
```

We now need to set up our new symmetric migration rates:

```{code-cell} python
m = 1e-3
set_migration_rates = [
    fwdpy11.SetMigrationRates(when=100, deme=0, migrates=[1 - m, m]),
    fwdpy11.SetMigrationRates(when=100, deme=1, migrates=[m, 1.0 - m]),
]
```

We now have the parts of our model.  To build a model, we need to create
an instance of {class}`fwdpy11.DiscreteDemography`:

```{code-cell} python
dmodel = fwdpy11.DiscreteDemography(
    mass_migrations=mass_migrations,
    migmatrix=migmatrix,
    set_migration_rates=set_migration_rates,
)
```

This class has a nice method for pretty-printing via `black`:

```{code-cell} python
print(dmodel.asblack())
```

Instances of {class}`fwdpy11.DiscreteDemography` are immutable, meaning that
the attributes are read-only and attempts to modify will raise exceptions:

```{code-cell} python
---
tags: [raises-exception]
---

dmodel.mass_migrations = None
```

However, you can get a {class}`dict` representation of the data, which you
can modify:

```{code-cell} python
import copy

dmodel_dict = copy.deepcopy(dmodel.asdict())
dmodel_dict["mass_migrations"] = None
dmodel2 = fwdpy11.DiscreteDemography(**dmodel_dict)
print(dmodel2.asblack())
```

:::{warning}

It is best practice to use {func}`copy.deepcopy` here. The
{class}`fwdpy11.DiscreteDemography` instances may contain
objects like {class}`numpy.ndarray` that only get copied by
reference when the `dict` is generated.  For safety/general
happiness, making a deep copy helps here.

:::

`fwdpy11` contains a small collection of pre-computed demographic models.
One of them is the common modeling scenario of two recently diverged populations.
Additionally, we supply some commonly-used models of human demography.
See {ref}`demographic-models` for details.

### A note of caution

:::{warning}

Programatically building up demographic models is tedious
and error-prone.  Errors end up in the literature, too (see
{cite}`Ragsdale2020-gl`).  Where possible, I **strongly** urge you
to test the correctness of a model against another implementation.
For example, one could compare results with no selection
and no recombination with the output of `msprime`. Another
approach is to compare to the predictions of
[moments](<https://bitbucket.org/simongravel/moments>) {cite}`Jouganous2017-tg`,
which is the approach taken by the `fwdpy11`
[statistical tests](<https://github.com/molpopgen/fwdpy11_statistical_tests>).

:::

Given the inherent difficulty of building up models, we hope to provide
a simpler approach in a future release.  There's a small group of people
currently hatching a plan to provide a common schema representing demographic
models of discrete demes.  The hope is that several pieces of software can
all use these schema.

(demographydebugger)=

### Debugging a demographic model

The parameters of a demographic model are checked at run time at two different places:

* Upon object construction.  The various event objects try to make sure that
  the parameter inputs are valid.
* If invalid events occur during a simulation,
  the simulation raises a `fwdpy11.DemographyError` exception.

It is clearly preferable for a simulation to detect errors as early as possible.
While bad inputs can be detected almost immediately, more subtle errors are only
detected during simulation, which may take a while.
A more efficient approach to checking your models is to use {class}`fwdpy11.DemographyDebugger`:

The class also generates a "report" with a verbal description of the model:

```{code-cell} python
d = fwdpy11.DemographyDebugger(
    initial_deme_sizes=[100],
    events=dmodel,
    simlen=150,
    deme_labels={0: "ANCESTRAL", 1: "DERIVED"},
)

print(d.report)
```

If the model were invalid, then an error would have been raised during initialization.

Let's take a look at passing in a bad model.  If we neglect to update
the migration rates, then the model considers there to be no ancestry specified for
deme `1`:

```{code-cell} python
---
tags: [raises-exception]
---
bad_dmodel = fwdpy11.DiscreteDemography(
    mass_migrations=mass_migrations, migmatrix=migmatrix
)
print(bad_dmodel.asblack())
d = fwdpy11.DemographyDebugger([100], bad_dmodel, 150, {0: "ANCESTRAL", 1: "DERIVED"})
```

Our `dmodel2` from above is also invalid, as deme `1`
never gets created since we deleted the mass migration event:

```{code-cell} python
---
tags: [raises-exception]
---
fwdpy11.DemographyDebugger([100], dmodel2)
```

(model-params)=

## Setting up the parameters for a simulation

Simulation parameters are stored in instances of {class}`fwdpy11.ModelParams`.
All of the examples shown above generated objects that we need to store in
such instances.

Like {class}`fwdpy11.DiscreteDemography`, {class}`fwdpy11.ModelParams` is
immutable after initialization.  Thus, I find it most convenient to
first put the objects into a {class}`dict` and then "explode" it to
initialize a {class}`fwdpy11.ModelParams`.

The following example uses what we discussed above:

```{code-cell} python
sregions = [fwdpy11.ExpS(beg=0.0, end=1.0, weight=1.0, mean=-0.2)]
recregions = [fwdpy11.PoissonInterval(beg=0.0, end=0.1, mean=1e-3)]
migmatrix = np.zeros(4).reshape(2, 2)
migmatrix[0, 0] = 1.0
mass_migrations = [
    fwdpy11.move_individuals(when=100, source=0, destination=1, fraction=0.5)
]
migrate = 1e-3
set_migration_rates = [
    fwdpy11.SetMigrationRates(when=100, deme=0, migrates=[1 - migrate, migrate]),
    fwdpy11.SetMigrationRates(when=100, deme=1, migrates=[migrate, 1.0 - migrate]),
]
dmodel = fwdpy11.DiscreteDemography(
    mass_migrations=mass_migrations,
    migmatrix=migmatrix,
    set_migration_rates=set_migration_rates,
)
dbg = fwdpy11.DemographyDebugger([100], dmodel, 150)
```

```{code-cell} python
gvalue = fwdpy11.Additive(
    scaling=2.0,
    gvalue_to_fitness=fwdpy11.GSSmo(
        [
            fwdpy11.Optimum(when=0, optimum=0.0, VS=1.0),
            fwdpy11.Optimum(when=100, optimum=1.0, VS=1.0),
        ]
    ),
)
pdict = {
    "nregions": [],
    "sregions": sregions,
    "recregions": recregions,
    "rates": (0.0, 1e-2, None),
    "gvalue": gvalue,
    "prune_selected": False,
    "simlen": 150,
    "demography": dmodel,
}
params = fwdpy11.ModelParams(**pdict)
```

There are a few new things here:

* The `nregions` field specifies where neutral mutations occur.  We leave it empty
  because we can add such mutations after the simulation is done.
* `prune_selected` tells the simulation what to do with fixed selected mutations
  after simplification. If `True`, they will be removed from the simulation.
  If `False`, they will be kept.  We keep them here because we are simulating an
  additive trait, and the fixed genetic background is part of the genetic value
  of an individual.  If we were instead simulating a standard population genetic model
  with multiplicative fitness effects, we could use `True` because all such models
  require is that relative fitnesses are preserved up to a multiplicative constant.
* The `rates` parameter is a list-like object with the neutral mutation rate,
  selected mutation rate, and the recombination rate, respectively.  See the
  next subsection for details.

Let's take a look at what we just built:

```{code-cell} python
print(params.asblack())
```

(model-params-rate-details)=

### Mutation and recombination rates

In general, the neutral mutation rate should be set to `0.0`.

The selected mutation rate is a non-negative {class}`float` representing the total mutation
rate per haploid genome per generation.

If you use instances of {class}`fwdpy11.Region` to set up the genetic map
(see {ref}`here <recregions>`), then you need to provide a non-negative {class}`float`
representing the total recombination rate per meiosis.  If you use instances
of {class}`fwdpy11.GeneticMapUnit` (see {ref}`here <geneticmapunit>`), then use
`None`.

(evolvets)=

## Running a simulation with tree sequence recording

To evolve the population with tree sequence recording, we make a call
to {func}`fwdpy11.evolvets`.  We also need an instance of {class}`fwdpy11.GSLrng`,
which is a random number generator.

```{code-cell} python
rng = fwdpy11.GSLrng(42)
pop = fwdpy11.DiploidPopulation(100, 1.0)
fwdpy11.evolvets(rng=rng, pop=pop, params=params, simplification_interval=25)
pop.generation
```

The `simplification_interval` parameter directs `fwdpy11` to run the tree
sequence simplification algorithm {cite}`Kelleher2018-mc` every 25 generations.  Simplifying
often uses less total memory but requires a longer total run time.  Simplifying too
infrequently may consume too much memory.  In general, I use a value of `100` and
try not to think about it too much.

:::{note}

This simulation involves one deme splitting into two. The way that
this model is written, the `gvalue` parameter applies to both
demes. See {ref}`here <mvdes>` for how to simulate mutations with
different effect sizes in different demes and {ref}`here <localadaptation>`
for a more complex example.

:::

:::{todo}
    How to track precise fixation times.
:::

(ancient-samples)=

## Recording "ancient samples" during a simulation

During a simulation, individuals can be marked as "preserved" or
"remembered".  This marking means that their nodes in the node
table and their metadata are maintained, allowing you to retrieve
the data for later analysis, reconstruct genotypes, etc..

To mark samples as preserved, you need to define a callable object
equivalent to the following function:

```{code-block} python

def record_samples(pop, sampler):
    pass

```

The first argument to this callable will be the instance of
{class}`fwdpy11.DiploidPopulation` that you've passed on to
{func}`fwdpy11.evolvets` (see {ref}`above <evolvets>`).
The second argument will be an instance of {class}`fwdpy11.SampleRecorder`.
This second argument will be created internally by `fwdpy11`
and passed to your callable.

The reason why we use the term "callable" here instead of "function"
is that you will probably want to write callable classes rather
than actual functions.  For example, imagine that we only want
to record individuals whose fitnesses are in the top `X%`
of those found in deme `Y` every `10th` generation:

```{code-cell} python
class GetTopFitnessess(object):
    def __init__(self, X, Y):
        self.X = X
        self.Y = Y

    def __call__(self, pop, sampler):
        if pop.generation % 10 == 0.0:
            md = np.array(pop.diploid_metadata, copy=False)
            inY = np.where(md["deme"] == self.Y)[0]
            w = md["w"][inY]
            q = np.quantile(w, self.Y)
            geqq = np.where(w >= q)[0]
            sampler.assign(inY[geqq])
```

The above class is not well-implemented at all, but it will suffice for our example,
where we pass an instance on to {func}`fwdpy11.evolvets`:

```{code-cell} python
rng = fwdpy11.GSLrng(42)
pop = fwdpy11.DiploidPopulation(100, 1.0)
recorder = GetTopFitnessess(0.5, 0)
fwdpy11.evolvets(
    rng=rng, pop=pop, params=params, simplification_interval=25, recorder=recorder
)
```

At the end of the simulation, we indeed have metadata from "ancient samples":

```{code-cell} python
amd = np.array(pop.ancient_sample_metadata, copy=False)
print(
    amd[
        -5:,
    ]
)
```

The following sampler types are built in to `fwdpy11`:

* {class}`fwdpy11.RandomAncientSamples`

(pop-from-tskit)=

## Tracking features of the population during a simulation

The concept of passing a callable into the simulation has more uses
than recording ancient samples.  Basically, we can track pretty
much anything.  The following example tracks the mean genetic value
in each deme every generation:

```{code-cell} python
class TrackGvalues(object):
    def __init__(self):
        self.data = []

    def __call__(self, pop, sampler):
        md = np.array(pop.diploid_metadata, copy=False)
        for i in [0, 1]:
            w = np.where(md["deme"] == i)[0]
            if len(w):
                self.data.append((pop.generation, i, md["g"][w].mean()))
```

```{code-cell} python
rng = fwdpy11.GSLrng(42)
pop = fwdpy11.DiploidPopulation(100, 1.0)
tracker = TrackGvalues()
fwdpy11.evolvets(
    rng=rng, pop=pop, params=params, simplification_interval=25, recorder=tracker
)
for i in tracker.data[-5:]:
    print(i)
```

(starting-from-msprime)=

## Starting a simulation with the output of a coalescent simulation

```{code-cell} python
import msprime

ts = msprime.simulate(200, Ne=100, recombination_rate=1e-3)
pop = fwdpy11.DiploidPopulation.create_from_tskit(ts)
pop.N
```

We can also start from a simulation of multiple demes:

```{code-cell} python
m = 1e-3
ts = msprime.simulate(
    Ne=100,
    recombination_rate=1e-3,
    population_configurations=[
        msprime.PopulationConfiguration(sample_size=100),
        msprime.PopulationConfiguration(sample_size=100),
    ],
    migration_matrix=np.array([0.0, m, m, 0.0]).reshape(2, 2),
)
pop = fwdpy11.DiploidPopulation.create_from_tskit(ts)
pop.N
md = np.array(pop.diploid_metadata, copy=False)
np.unique(md["deme"], return_counts=True)
```

See {ref}`here <precapitation>` for a detailed example.

Currently, `fwdpy11` does **not** support starting from `msprime` simulations
with mutations.  The variants will simply be ignored.

:::{note}

When simulating very large genomic regions, it is not a good
idea to use the standard model in `msprime`. See {cite}`Nelson2020-my`
for details and also {class}`msprime.DiscreteTimeWrightFisher`.
Also note that starting with `msprime` output doesn't
guarantee that you don't need to "burn in" the simulation
(see {ref}`here <howlongtorun>`).

:::

:::{warning}

Not all genetic maps supported by `fwdpy11` are possible
in `msprime`. Specifically, discrete "jumps" in the
genetic map (*e.g.*, {class}`fwdpy11.BinomialPoint`) are
not compatible with `msprime`.  Further, the approach
outlined [here](<https://msprime.readthedocs.io/en/stable/tutorial.html#multiple-chromosomes>)
is not sufficient to get the recombination rates right
between separate "chromosomes".

:::

(tskittransfer)=

## Converting the data to a {class}`tskit.TreeSequence`

The tree sequence data structures may be converted to the analogous `tskit`
objects using {func}`fwdpy11.DiploidPopulation.dump_tables_to_tskit`,
which returns a {class}`tskit.TreeSequence`.

The most basic usage is:

```{code-block} python

ts = pop.dump_tables_to_tskit()

```

When you have the data stored as a {class}`tskit.TreeSequence`,
information about individuals, mutations, etc., is stored as table metadata.
See {ref}`here <tskit_metadata_vignette>` to learn how to decode the metadata.

:::{note}

Once `tskit` 0.3.0 is released, the metadata encoding will change
quite a bit and it will be simpler and more efficient to decode.

:::

You may provide a `dict` that reflects the simulation parameters used.  This `dict`
will be part of the provenance information encoded in the {class}`tskit.TreeSequence`.
For example:

```{code-block} python

# Assuming mp is a fwdpy11.ModelParams:
ts = pop.dump_tables_to_tskit(parameters={"model": str(mp), "seed": 12345})

```

Ultimately, it is up to you to decide what to include in `parameters`.
For example, it could be a script:

```{code-block} python

# Bonus points for somehow including the git commit hash corresponding
# to the version of the script that you used !
parameters = {"script": "/path/to/script", "type": "script", "seed": 1234}

```

You can go further than that, and even include the entire script.
It turns out that Python files know who they are and can read themselves:

```{code-block} python

def read_self():
    with open(__file__, "r") as f:
        script = f.read()
    return script


script = read_self()
parameters = {"script": script, "type": "script", "seed": 1234}

```

In order to get the provenance information back out from a {class}`tskit.TreeSequence`:

```{code-block} python

import json

provenance = json.loads(ts.provenance(0).record)

```

If you recorded an instance of {class}`fwdpy11.ModelParams` as your `parameters`, you
can even reconstruct the original object (if you have the correct modules imported).
For example, if we assume that we encoded the model parameters as shown two listings ago:

```{code-block} python

import tskit
import json
import numpy as np
import fwdpy11

ts = tskit.load("sim.trees")
provenance = json.loads(ts.provenance(0).record)
array = np.array  # Annoyance!
tmp = eval(provenance["parameters"]["model"])

```

It is possible for model parameters to contain `numpy` arrays.  Unfortunately, their
string representations are not namespace-qualified, meaning that they say `array` rather
than `numpy.array` or `np.array`. Thus, I made a type alias so that the `eval` would work.

(savingsimstodisk)=

## Saving the results of a simulation to disk

To dump the output of a simulation to an uncompressed binary file, use
{func}`fwdpy11.DiploidPopulation.dump_to_file`:

```{code-block} python

pop.dump_to_file("pop.bin")

```

To restore a population from such a file, call the static method
{func}`fwdpy11.DiploidPopulation.load_from_file`:

```{code-block} python

pop = fwdpy11.DiploidPopulation.load_from_file("pop.bin")

```

There are two means of storing a population in {mod}`pickle` format.
The first is:

```{code-block} python

import gzip
import pickle

with gzip.open("pickled_pop.gz", "wb") as f:
    pickle.dump(pop, f)

```

However, this method is not the most efficient.  The following takes less
memory, which is probably important for really big simulations:

```{code-block} python

with gzip.open("pickled_pop.gz", "wb") as f:
    pop.pickle_to_file(f)

```

:::{warning}

These two pickling methods create very different files!
The first can be read back in with {func}`pickle.load`
but the second **must** be read back in with the
static function {func}`fwdpy11.DiploidPopulation.load_from_pickle_file`.

:::

### Outputting a `tskit` "trees file"

```{code-block} python

ts = pop.dump_tables_to_tskit()
ts.dump("treefile.trees")

```


