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
(softselection)=

# Soft selection with discrete demes

This page describes implementing models of demographic events affecting
multiple demes.  This functionality was first released in version 0.6.0
and makes use of low-level types added in 0.5.3.

## Overview

This chapter describes how to generate detailed demographic models of individuals
evolving withing discrete sub-populations, or demes.

If you want to skip the details and see basic models provided by `fwdpy11`, then
see {ref}`demographic-models`.

Once you've digested this section, you may want to read {ref}`mvdes`.

### The model

The model here is soft selection {cite}`Levene1953-yn`, meaning that the number of
breeding individuals ("adults") in each deme is fixed at a certain value.
A nice overview of this model and how it compares to others in {cite}`Felsenstein1976-sb`.
You may also find {cite}`Christiansen1974-ac` and {cite}`Christiansen1975-ch` useful.

Each generation, we generate offspring ("juveniles") in each deme.  In the absence of
migration, all parents of all offspring come from the offspring deme.  With a migration
matrix, we first choose the "source" deme, pick parents from that deme, and then
create the offspring in the offspring deme.  Thus, we are modeling juvenile migration.
Selfing is a property of demes and is applied to parents: we choose a parental deme,
then choose a parent, and then decide if that parent selfs.

:::{note}

The migration behavior changed in 0.6.2!  (This is mainly a note for
the developers.)

:::

### The timings of events

Below, we discuss various events that may happen.  These event types
include things like deme size changes, "mass migration" events, etc..
These events will occur at a certain time in a simulation. That time
refers to the birth time of a generation and the events are applied
*prior* to generating offspring, meaning that the events happen *to
the parents*.  For example, if half of deme zero moves and colonizes
a new deme (deme 1), then that means that half of the current alive individuals
(possible parents) have their `deme` field changed from zero to one
prior to generating any offspring.

The objects events all have an attribute called `when`, which
parameterizes when the event occurs in a simulation.  The value of `when`
is with respect to the current generation time of the population
({attr}`fwdpy11.DiploidPopulation.generation`). When events are registered
to occur prior to this time, you will see warnings.  See {func}`fwdpy11.evolvets`
for details.

Events scheduled for prior than the population's current time are allowed
because there are valid modeling reasons to allow them.  For example, you
may want to evolve for a while and then change some other model parameter
like the recombination rate, and then keep evolving.

(soft-sel-deme-setup)=

## Setting the initial demes in a simulation

At the start of a simulation, you may assign diploids to demes
when constructing an instance of {class}`fwdpy11.DiploidPopulation`.
For example, to initialize a population with 25 individuals in demes `0` and `1`:

```{code-cell} python
import fwdpy11
import numpy as np

pop = fwdpy11.DiploidPopulation([25, 25], 1.0)
md = np.array(pop.diploid_metadata, copy=False)
pop.deme_sizes()
for m in pop.diploid_metadata:
    for n in m.nodes:
        assert m.deme == pop.tables.nodes[n].deme
```

Another method involves mass migration events at the beginning of a simulation.
See {ref}`massmigrations`.

## The DiscreteDemography class

The demographic events are stored in instances of {class}`fwdpy11.DiscreteDemography`.
These events, whose interface is described below, are passed in `list` objects
when creating a {class}`fwdpy11.DiscreteDemography` instance.

These instances may be used to parameterize the `demography` field of a
{class}`fwdpy11.ModelParams` instance.  To illustrate this, here is a
function that we'll use repeatedly below:

```{code-cell} python
def setup_and_run_model(pop, ddemog, simlen, recorder=None, seed=654321):
    pdict = {
        "nregions": [],
        "sregions": [],
        "recregions": [],
        "rates": (
            0,
            0,
            0,
        ),
        "gvalue": fwdpy11.Multiplicative(2.0),
        "demography": ddemog,
        "simlen": simlen,
    }
    params = fwdpy11.ModelParams(**pdict)
    rng = fwdpy11.GSLrng(seed)
    fwdpy11.evolvets(rng, pop, params, 100, recorder)
```

We will also define a simple class to record all deme sizes over time:

```{code-cell} python
class SizeTracker(object):
    def __init__(self):
        self.data = []

    def __call__(self, pop, sampler):
        self.data.append((pop.generation, pop.N, pop.deme_sizes()))
```

### Compatibility with previous versions of fwdpy11

:::{versionchanged} 0.8.0

You now must specify `simlen` manually.

:::

Previous versions only supported size changes within a single deme.  These size changes were
parameterized via a `numpy` array specifying the size at each time point.  It is still possible
to specify the demography using that approach:

```{code-cell} python
N = np.array([10] * 10 + [5] * 5 + [10] * 10, dtype=np.uint32)
pdict = {
   "nregions": [],
   "sregions": [],
   "recregions": [],
   "rates": (
       0,
       0,
       0,
   ),
   "gvalue": fwdpy11.Multiplicative(2.0),
   "demography": fwdpy11.DiscreteDemography(set_deme_sizes=N),
   "simlen": len(N),
}
params = fwdpy11.ModelParams(**pdict)
rng = fwdpy11.GSLrng(654321)
pop = fwdpy11.DiploidPopulation(10, 1.0)
fwdpy11.evolvets(rng, pop, params, 100)
```

## Event types

The following sub-sections describe the various types of demographic
events allowed during a simulation.

(massmigrations)=

### Mass migrations

Mass migration events represent the "bulk" movement of individuals
in a single generation.  Such events allow you to model population
splits, merges, etc..

These events are represented by instances
of {class}`fwdpy11.MassMigration`.  Currently, you create instances
of this type using one of the following two functions:

* {func}`fwdpy11.copy_individuals`
* {func}`fwdpy11.move_individuals`

As the name implies, the first function creates an event that *copies*
individuals from a source deme to a destination.  The latter *moves*
them.

Both functions take five arguments, which may be used either named
or unnamed.  In order, they are:

* `when`: the time (generation) when the event will occur
* `source`: the ID of the source deme
* `destination`: the ID of the destination deme
* `fraction`: the fraction (proportion) of `source` moved/copied to `dest`.
* `resets_growth_rate`: If `True`, the event resets the growth rate to {attr}`fwdpy11.NOGROWTH`
  in **both** `source` and `dest`. If `False`, growth rates remain unchanged.
  The default is `False`.

:::{note}

When a mass migration event *copies* individuals from deme,
the individuals copied are sampled *without replacement*.  Thus,
if the fraction copied is 1.0, then every individual is copied.

:::

These operations act on proportions of populations rather than on numbers
of individuals. Multiple events in a single generation are allowed, see
{ref}`multiple-mass-migrations`.

#### Setting the initial state of a simulation

Let's look at an example where we use mass migration events to set up
"who is where" at the start of a simulation.  Since events happen in
the *parental* generation, we can use mass migrations to set up
what demes individuals are in by applying events at generation 0.

The main difference between this method and that shown in
{ref}`soft-sel-deme-setup` is that these events move or copy *random*
individuals to new demes whereas using the  `__init__` approach
builds the individuals in each deme sequentially.

For example, if we wish to start a simulation with 50 individuals in
demes 0 and 50 in deme 1, we have two options:

1. Start with 50 individuals and *copy* them to deme 1 in generation 0
2. Start with 100 individuals and *move half of* them to deme 1 in generation 0

Here is the version implemented via a  copy:

```{code-cell} python
pop = fwdpy11.DiploidPopulation(50, 1.0)
copy = [fwdpy11.copy_individuals(when=0, source=0, destination=1, fraction=1.0)]
ddemog = fwdpy11.DiscreteDemography(mass_migrations=copy)
setup_and_run_model(pop, ddemog, 1)
pop.deme_sizes()
```

Here is what our object looks like:

```{code-cell} python
copy[0]
```

Here is the version using a move:

```{code-cell} python
pop = fwdpy11.DiploidPopulation(100, 1.0)
move = [fwdpy11.move_individuals(0, 0, 1, 0.5)]
ddemog = fwdpy11.DiscreteDemography(mass_migrations=move)
setup_and_run_model(pop, ddemog, 1)
pop.deme_sizes()
```

For comparison, here is the object specifying the move:

```{code-cell} python
move[0]
```

(multiple-mass-migrations)=

#### Multiple mass migrations

To specify multiple events, simply add more events to your list.
The events to not have to be sorted in any specific way.  Any sorting
requirements get handled internally.

Multiple events involving the same source population in the same generation
need some explaining.   If the events are copies, things will tend to "just
work":

```{code-cell} python
pop = fwdpy11.DiploidPopulation(50, 1.0)
copy = [fwdpy11.copy_individuals(0, 0, 1, 1.0), fwdpy11.copy_individuals(0, 0, 2, 1.0)]
ddemog = fwdpy11.DiscreteDemography(mass_migrations=copy)
setup_and_run_model(pop, ddemog, 1)
pop.deme_sizes()
```

When the events are moves, it is not possible to move more than 100%
of the individuals.  Attempting to do so will raise a `ValueError`
exception:

```{code-cell} python
pop = fwdpy11.DiploidPopulation(50, 1.0)
# Move all of deme 0 into demes 1 and 2,
# which means we're trying to move 200%
# of deme 0...
move = [fwdpy11.move_individuals(0, 0, 1, 1.0), fwdpy11.move_individuals(0, 0, 2, 1.0)]
# ... which is not allowed
try:
    ddemog = fwdpy11.DiscreteDemography(mass_migrations=move)
except ValueError as e:
    print(e)
```

#### The rate of drift

Moving versus copying individuals is an important modeling choice.
When you move individuals from one deme to another, the rate of drift
changes in the source deme (as its size is reduced).  This reduction
in size is also a sudden bottleneck.

Copying, on the other hand, does not change the rate of drift in the source
deme.  However, it does seem to imply some sudden increase in fecundity that
both came from nowhere and was short-lived.

(set-deme-sizes)=

### Instantaneous deme size changes

Instantaneous changes in deme size are managed by instances of
{class}`fwdpy11.SetDemeSize`.

This class is relatively straightforward to use, so let's dive right in:

```{code-cell} python
pop = fwdpy11.DiploidPopulation([20, 20], 1.0)
dd = fwdpy11.DiscreteDemography(
    set_deme_sizes=[fwdpy11.SetDemeSize(when=5, deme=1, new_size=100)]
)
st = SizeTracker()
setup_and_run_model(pop, dd, 10, st)
for i in st.data:
    print(i)
```

You may also kill off demes by setting their size to zero:

```{code-cell} python
pop = fwdpy11.DiploidPopulation([20, 20, 20], 1.0)
dd = fwdpy11.DiscreteDemography(
    set_deme_sizes=[fwdpy11.SetDemeSize(when=5, deme=1, new_size=0)]
)
st = SizeTracker()
setup_and_run_model(pop, dd, 6, st)
for i in st.data:
    print(i)
```

### Changing growth rates

Instances of {class}`fwdpy11.SetExponentialGrowth` manage the exponential growth rates per deme.
Growth rates less than one indicate population decline, greater than one means growth
and {attr}`fwdpy11.NOGROWTH` is equal to 1.0 to indicate no growth.

Let's look at an example:

```{code-cell} python
pop = fwdpy11.DiploidPopulation([50], 1.0)
g = [fwdpy11.SetExponentialGrowth(when=0, deme=0, G=1.1)]
dd = fwdpy11.DiscreteDemography(set_growth_rates=g)
st = SizeTracker()
setup_and_run_model(pop, dd, 6, st)
for i in st.data:
    print(i)
```

The deme sizes each generation must be integer values.  The simulation uses C/C++ rules for
rounding double-precision values to integer values. The function `numpy.rint` uses the same
rules:

```{code-cell} python
N0 = np.float64(50.0)
for i in range(6):
   Ni = N0 * np.power(1.1, i + 1)
   print(i + 1, Ni, np.rint(Ni))
```

You may need to keep the rounding policy in mind when trying to predict final deme sizes when testing
or when trying to convert a model from continuous time into discrete time.

### Changing the selfing rate

Instances of {class}`fwdpy11.SetSelfingRate` affect the rate of selfing-versus-outcrossing in different
demes, or to change the rate within a deme over time. The default is that individuals don't self
unless they are picked twice as a parent by chance.

Using this type is straightforward.  Before we dive in, we will create a new recorder
type to track parents each generation:

```{code-cell} python
class ParentTracker(object):
    def __init__(self):
        self.data = []

    def __call__(self, pop, sampler):
        for i in pop.diploid_metadata:
            self.data.append((i.label, i.deme, i.parents))
```

Let's run a simulation for a couple of generations:

```{code-cell} python
pop = fwdpy11.DiploidPopulation([5, 5], 1.0)
sr = [fwdpy11.SetSelfingRate(when=0, deme=1, S=1.0)]  # Deme 1 always selfs
dd = fwdpy11.DiscreteDemography(set_selfing_rates=sr)
pt = ParentTracker()
setup_and_run_model(pop, dd, 2, pt)
```

In our output, the deme label is the second value in each tuple, and any individual
in deme 1 has the same parent listed twice because they were the product of a selfing event:

```{code-cell} python
for i in pt.data:
    print(i)
```

(In the above output, the parent IDs are the indexes of the parental individuals from their
generation.)

(migration)=

### Migration

For models with multiple demes, migration between then is managed by an
instance of {class}`fwdpy11.MigrationMatrix`.

For a migration matrix `M`, the default interpretation of `M[i, j]` is the
fraction of deme `i` that will be replaced by migrations from deme `j`. The
entry `M[i, i]` represents the non-migrant fraction of deme `i`'s ancestry.
The matrix is "row-major" meaning that rows refer to migration into source demes.
This definition of the migration matrix corresponds to that found in several
different sources ({cite}`Christiansen1974-ac`, {cite}`Christiansen1975-ch`).
This definition of migration is also what diffusion models assume (*e.g.* {cite}`Jouganous2017-tg`)
as well as coalescent simulations like *msprime* {cite}`Kelleher2016-cb`.

For example, consider the following matrix:

```{code-cell} python
m = np.array([0.9, 0.1, 0.5, 0.5]).reshape(2, 2)
m
```

The first row corresponds to the ancestry of deme `0`, such that 90% of offspring will be
non-migrants and 10% will be migrants from deme `1`:

```{code-cell} python
m[
   0,
]
```

To be concrete, if the size of deme `0` in the next generation is 1,000, then the expected
number of migrant and non-migrant offspring of offspring in deme `0` is:

```{code-cell} python
m[
   0,
] * 1e3
```

The second row implies that half the ancestry of deme `1` is due to migrants and half
due to non-migrants:

```{code-cell} python
m[
   1,
]
```

The `numpy` array is sufficient to construct our demographic model:

```{code-cell} python
d = fwdpy11.DiscreteDemography(migmatrix=m)
d.migmatrix
```

By default, there is no migration, which is represented by the value `None`.  For example,
the following model has no migration events:

```{code-cell} python
# Define demographic events w/o any migration stuff
d = fwdpy11.DiscreteDemography(set_deme_sizes=[fwdpy11.SetDemeSize(0, 1, 500)])
d.migmatrix is None
```

In order to specify a model with no initial migration, you may use an identity matrix:

```{code-cell} python
d = fwdpy11.DiscreteDemography(migmatrix=np.identity(2))
d.migmatrix
```

The only reason to use the identity matrix is to start a simulation with no migration
and then change the rates later via instances of {class}`fwdpy11.SetMigrationRates`.
To see this in action, we'll first generate a new type to track if parents of
offspring in deme 1 are migrants or not:

```{code-cell} python
class MigrationTracker(object):
    def __init__(self, N0):
        self.N0 = N0
        self.data = []

    def __call__(self, pop, sampler):
        for i in pop.diploid_metadata:
            if i.deme == 1:
                p = []
                for j in i.parents:
                    if j < self.N0:
                        p.append((j, True))
                    else:
                        p.append((j, False))
                self.data.append((pop.generation, i.label, i.deme, p))
```

```{code-cell} python
# No migration at first
mm = np.identity(2)
# In generation 3, reset migration rates for deme 1 such
# that parents are equally likey from both demes.
cm = [fwdpy11.SetMigrationRates(3, 1, [0.5, 0.5])]
dd = fwdpy11.DiscreteDemography(migmatrix=mm, set_migration_rates=cm)
pop = fwdpy11.DiploidPopulation([10, 10], 1.0)
mt = MigrationTracker(10)
setup_and_run_model(pop, dd, 4, mt)
```

```{code-cell} python
for i in mt.data:
    nmig = 0
    if i[1] > 10:
        if i[3][0][1] is True:
                nmig += 1
        if i[3][1][1] is True:
            nmig += 1
    mstring = ""
    if nmig > 0:
        mstring = "<-  migrant parent".format(nmig)
    if nmig > 1:
        mstring += "s"
        print(i, mstring)
```

#### An alternative model of migration

The description of migration rates above implies that migration events are
independent of of source deme sizes.  To revisit our earlier example:

```{code-cell} python
m = np.array([0.9, 0.1, 0.5, 0.5]).reshape(2, 2)
# The is the expected number of parents from demes 0 and 1
# to offspring born in deme 0:
m[
   0,
] * 1000
```

`fwdpy11` allows for a different migration scheme where the size of the source deme
matters.  For this model, `M[i ,j]` is the probability that an individual with parents from
deme `j` is born in deme `i`.  Internally, the migration matrix entries
`M[i, j]` are multiplied by the size of the *source* demes, which means that
larger demes with nonzero migration rates to other demes have a larger chance
of being sources of migrant offspring.

For example:

```{code-cell} python
deme_sizes = np.array([1000, 2000])
m
md = m * deme_sizes
# The following line divides each
# row by its sum
md / np.sum(md, axis=1)[:, None]
```

The first matrix is the same as in the preceding section--90% of the offspring in deme
`0` will have parents from deme `0`.  In the second matrix, that fraction is reduced to
about 82% because deme `1` is twice as large as deme `0`.

To enable this migration model, create an instance of {class}`fwdpy11.MigrationMatrix` and
pass `True` as the second parameter:

```{code-cell} python
M = fwdpy11.MigrationMatrix(m, True)
d = fwdpy11.DiscreteDemography(migmatrix=M)
```

This will also work, but is less explicit:

```{code-cell} python
d = fwdpy11.DiscreteDemography(migmatrix=(m, True))
```

:::{note}

This model of migration will typically give *different* results
from diffusion models and coalescent simulations!

:::

(migration-and-selfing)=

## Examples of models

### Isolation with migration, or "IM"

Consider two demes that split apart `T` time units ago and then grow to different
sizes in the present.  After the split, migration occurs between the two demes. The
demographic model has the following parameters:

* `Nanc`, the ancestral population size.
* `T`, the time of the split, which is in units of `Nanc`.
* `psplit`, the proportion of the ancestral population that splits off to found deme `1`.
* `N0`, the final size of deme `0`, relative to `Nanc`.
* `N1`, the final size of deme `1`, relative to `Nanc`.
* `m01`, the migration rate from deme `0` to deme `1`.
* `m10`, the migration rate from deme `1` to deme `0`.

Here is the model in its entirety, with no mutation and no recombination.
First, we will set up the demographic events.  The population with evolve
for `Nanc` generations before the split.

```{code-cell} python
Nanc = 100
T = 0.2
psplit = 0.33
N0, N1 = 2, 3
m01, m10 = 0.01, 0.0267

# The split event
split = [fwdpy11.move_individuals(when=Nanc, source=0, destination=1, fraction=psplit)]
# Get growth rates and set growth rate changes,
# taking care to handle our rounding!
gens_post_split = np.rint(Nanc * T).astype(int)
N0split = np.rint(Nanc * (1.0 - psplit))
N0final = np.rint(N0 * Nanc)
N1split = np.rint(Nanc * psplit)
N1final = np.rint(N1 * Nanc)
G0 = fwdpy11.exponential_growth_rate(N0split, N0final, gens_post_split)
G1 = fwdpy11.exponential_growth_rate(N1split, N1final, gens_post_split)
growth = [
    fwdpy11.SetExponentialGrowth(Nanc, 0, G0),
    fwdpy11.SetExponentialGrowth(Nanc, 1, G1),
]

# Set up the migration matrix for two demes, but only
# deme zero exists.
m = fwdpy11.migration_matrix_single_extant_deme(2, 0)
# The rows of the matrix change at the split:
cm = [
    fwdpy11.SetMigrationRates(Nanc, 0, [1.0 - m10, m10]),
    fwdpy11.SetMigrationRates(Nanc, 1, [m01, 1.0 - m01]),
]
d = fwdpy11.DiscreteDemography(
    mass_migrations=split, set_growth_rates=growth, set_migration_rates=cm, migmatrix=m
)
```

The above code made use of two helper functions:

* {func}`fwdpy11.exponential_growth_rate`
* {func}`fwdpy11.migration_matrix_single_extant_deme`

Finally, we can run it:

```{code-cell} python
pop = fwdpy11.DiploidPopulation(Nanc, 1.0)
setup_and_run_model(pop, d, Nanc + gens_post_split)
```

Now we check the final population sizes and make sure they are correct:

```{code-cell} python
ds = pop.deme_sizes()
assert ds[1][0] == N0final
assert ds[1][1] == N1final
```

This model is common enough that you shouldn't have to implement it from
scratch each time.  For this reason, we provide it in {func}`fwdpy11.demographic_models.IM.two_deme_IM`.

```{code-cell} python
import fwdpy11.demographic_models.IM

dmodel = fwdpy11.demographic_models.IM.two_deme_IM(
    Nanc, T, psplit, (N0, N1), (m01, m10), burnin=1.0
)
pop2 = fwdpy11.DiploidPopulation(Nanc, 1.0)
setup_and_run_model(pop2, dmodel.model, dmodel.metadata.simlen)
assert pop.generation == pop2.generation
assert pop2.generation == dmodel.metadata.simlen
ds2 = pop2.deme_sizes()
assert np.array_equal(ds[0], ds2[0])
assert np.array_equal(ds[1], ds2[1])
```

See {ref}`IMexample` for an example of using this function to compare to results
from diffusion models.


