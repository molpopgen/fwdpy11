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

(advancedtopics)=

# Advanced topics

## Executing multiple replicates in parallel using `concurrent.futures`

You have many options when it comes to running many replicates of a model.
For example, you could write a single Python script plus a shell script
that your local cluster understands, sending your simulation to hundreds
of compute nodes.

This section introduces {mod}`concurrent.futures` as a method to execute
many replicates in separate Python *processes* using a single script.

Fully worked-out examples of varying complexity can be found here:

* {ref}`bgs_vignette`
* {ref}`IMexample`
* {ref}`migtest`
* {ref}`precapitation`
* {ref}`recapitation`

A common set of idioms used by these examples are:

* The function to run the simulation typically does not return
  the instance of {class}`fwdpy11.DiploidPopulation`.  Such
  a return requires pickling, which is very expensive.  If you
  need to store these population, write them to files rather
  than returning them. See {ref}`here <savingsimstodisk>`.
* If you do return something to the main "collector"
  process, keep it simple!
* Random number seeds are set in `__main__` using
  {func}`numpy.random.randint` after calling
  {func}`numpy.random.seed` using a user-provided seed.

(howlongtorun)=

## How long to run the simulation

Deciding how long to run a simulation depends on whether or not
a "burn in" phase is required to get the model to statistical
equilibrium.  For example, if you want to get the steady-state
properties of a model with many selected mutations and linkage,
you will need to burn in.  If you want to know how the short term
dynamics of the neutral site-frequency spectrum are affected by suddenly
introducing strongly-beneficial mutations, then you may not need to
burn in.

So, then, how long to burn in?  Clearly, the answer should be "the
least that you have to"!  As always, there are a few different
things to consider:

* You want the *distribution* of final outcomes to be "correct".
  Theory often gives us expressions for expectations and
  sometimes for variances of things we care about, like
  the number of mutations in a sample. However, comparing
  distributions will often require a statistical analysis
  comparing the output of independent simulations.
  For example, for a neutral model, the distribution of the
  number of mutations under an infinitely-many sites model
  for a small sample should match `msprime` very well.
  For a steady-state model of recurrent hitch-hiking
  {cite}`Kaplan1989-rt`, results from small samples
  should match coalescent simulations such as {cite}`Kern2016-wg`.
  Clearly, all of the different methods
  should match analytical results where possible.
* You probably want all of your trees to be fully coalesced
  to a single common ancestor.  The time back to a final
  ancestor has a very long tail.  In models with recombination,
  it is not uncommon to have a tree or three with multiple roots
  at the end of a simulation.

For the first point, a burn-in of `10N` generations or so
is usually sufficient. Before we simulated with tree sequence
recording, we simulated entire "genomes" (neutral mutations and all),
which was rather slow.  (Here "we" is the field in general.)  We then
(hopefully!) compared our results to something like `msprime`.  The
agreement was really good, so one answer is "about `10N`.  I'm being
vague about `N` here--typically, it would be equal to the starting effective
population size of your population (*e.g.* in the absence of selection).

With tree sequence recording, we run into the "uncoalesced marginal
trees" issue mentioned above.  We could just simulate much longer
to get rid of this problem, but it is simply easier to start
with a tree sequence from `msprime` (see {ref}`here <starting-from-msprime>`).

(eyre-walker)=

## The "Eyre-Walker" model of complex traits

{cite}`Eyre-Walker2010-rs` describes a model relating mutations with
fitness effect {math}`S = |2Ns|` to {math}`z`, their effect on
a phenotype/trait according to {math}`z = \delta S^\tau(1+\epsilon)`.
Here, {math}`\epsilon` is a draw from a Gaussian distribution with mean zero,
{math}`\delta` is {math}`1` or {math}`-1` with equal probability, and
{math}`\tau` "tunes" the correlation between fitness effect and effect on the trait.

Implementing this model is quite straightforward, as the trait values do not affect
the dynamics of the model.  For this flavor of a "pure pleiotropy" model,
the technical details reduce to a standard population genetic model
where we tack on the trait values at the end.

Let's set up a model where {math}`S` is exponentially distributed with mean `-20`.
We'll run a small population for a few generations.  This model will be nowhere
near equilibrium, but we're just using it as an example:

```{code-cell} python
import fwdpy11
import numpy as np

N = 1000
sregions = [fwdpy11.ExpS(beg=0.0, end=1.0, weight=1.0, mean=-20, scaling=2 * N)]
recregions = [fwdpy11.PoissonInterval(beg=0.0, end=0.1, mean=1e-3)]
gvalue = fwdpy11.Multiplicative(scaling=2.0)
pdict = {
    "nregions": [],
    "sregions": sregions,
    "recregions": recregions,
    "rates": (0.0, 1e-3, None),
    "gvalue": gvalue,
    "prune_selected": False,
    "simlen": 150,
}
params = fwdpy11.ModelParams(**pdict)
pop = fwdpy11.DiploidPopulation(N, 1.0)
rng = fwdpy11.GSLrng(54321)
fwdpy11.evolvets(rng, pop, params, 100)
print(len(pop.tables.mutations))
```

Define a function relating {math}`S` to {math}`z`:

```{code-cell} python
def getz(S, tau, sigma):
    if np.random.uniform() < 0.5:
        delta = 1
    else:
        delta = -1
    epsilon = np.random.normal(loc=0, scale=sigma, size=1)[0]
    return delta * np.power(np.abs(S), tau) * (1.0 + epsilon)
```

Apply the function and look at the results:

```{code-cell} python
np.random.seed(101010)
zvals = {}
for i, m in enumerate(pop.tables.mutations):
    zvals[m.key] = getz(2 * N * pop.mutations[m.key].s, 0.5, 1.0)

for k, v in zvals.items():
    print(pop.mutations[k].s, v)
```

Traverse the tree sequence to get individual phenotypes under a
strictly additive model:

```{code-cell} python
phenotypes = np.zeros(pop.N)
node_to_individual = {}
for i, j in enumerate(pop.diploid_metadata):
    assert j.nodes[0] not in node_to_individual
    assert j.nodes[1] not in node_to_individual
    node_to_individual[j.nodes[0]] = i
    node_to_individual[j.nodes[1]] = i
ti = fwdpy11.TreeIterator(pop.tables, pop.alive_nodes, update_samples=True)
for t in ti:
    for m in t.mutations():
        for n in t.samples_below(m.node):
            phenotypes[node_to_individual[n]] += zvals[m.key]
```

The trait value distribution is:

```{code-cell} python
np.unique(phenotypes, return_counts=True)
```

The mean trait value and the genetic variance are:

```{code-cell} python
phenotypes.mean()
phenotypes.var()
```

For our final trick, let's store them in the individual metadata:

```{code-cell} python
md = np.array(pop.diploid_metadata, copy=False)
md["g"][:] = phenotypes

for i, j in zip(md["g"][:5], phenotypes[:5]):
    print(i, j)
```

The trick is that the `numpy` array is a *non-owning* array, meaning
that it is a simple "view" of the underlying `C++` data.  Further,
it happens to be read/write, and thus we can modify it.  (Try not
to abuse this.  More often than not, you'll just break stuff.)

I feel that some comments about this model are warranted:

* The model is a simplification of earlier work by Keightley and Hill
  ({cite}`Keightley1988-rn`, {cite}`Keightley1990-az`), which should be cited
  alongside {cite}`Eyre-Walker2010-rs`.
* {cite}`Keightley1988-rn` and {cite}`Keightley1990-az`) show that this model predicts heritability
  (the genetic variance) is linear-ish with `N`, with the details depending
  somewhat on the model parameters.  The biological reasonableness
  of that prediction is dubious.  {cite}`Johnson2005-bp` discuss this point
  in some detail, and give other relevant references.
* Depending somewhat on the parameters of the `getz` function, one will eventually
  generate an intermediate-frequency variant with a massive {math}`z`,
  meaning that it would explain a considerable amount of the genic
  variance for the trait (high {math}`2pqz^2`).  Such outcomes are contrary
  to the results of human GWAS.
* It becomes tempting to start doing arbitrary things when applying this model
  to simulation output.  For example, only treating mutations in certain frequency
  ranges as affecting trait values.  However, doing so means the predictions
  made by such procedures are not natural outcomes of evolutionary models.
  Consider only applying the `getz` function to mutations with frequency
  {math}`< x` in order to say something about "rare alleles".
  This treatment of the data actually changes the model: mutations that are common
  now (at the end of the simulation) were rare at some point in the past,
  and had frequencies {math}`< x`.  Therefore, the model is one where mutations
  suddenly stop affecting the trait once they hit some critical frequency.
  Presumably, the same variants would affect the trait again should they drift to
  a frequency below {math}`x` if the simulation is run a bit longer.

(gvalues_python)=

## Writing genetic value classes in Python

Version `0.9.0` added the ability to write custom genetic value types in
Python.  Such classes are useful for prototyping new models
and/or generating example models for teaching/illustrating concepts.

:::{versionadded} 0.9.0

:::

### Concepts

(python_genetic_value_dimensions)=

#### The "dimensionality" of a model

`fwdpy11` supports both univariate and multivariate trait simulations.
Thus, a genetic value implicitly has a dimensionality. For the
case of a standard population-genetic simulation where mutations
directly affect fitness, that dimensionality is `1`.  Simply put,
the genetic value is fitness.

For a simulation of mutations affecting `n` traits, the dimensionality
is `n`.

The dimensionality matters because we need to define a single genetic value
to populate the {attr}`fwdpy11.DiploidMetadata.g` field for an offspring.
Doing so is simple for a one-dimensional trait.  For a simulation with
pleiotropy, we have a decision to make.  Should {attr}`fwdpy11.DiploidMetadata.g`
store the Euclidean distance from the optima?  Or should it store the
individual's value for some "focal" trait?  The `API` described here is
flexible enough to accommodate both.

The dimensionality of your genetic value to fitness map and your genetic value
object must match.  For the trait examples discussed above, this requirement
should make good sense.  For the case of modeling correlated effect sizes
of mutations in different demes (see {ref}`here <mvdes>`), we still require
this although it makes less sense. (Such simulations are limited to either
the standard population genetic case or the case of selection on a single
trait.  But relaxing this requirement for some models but not others becomes
problematic.)

The current version of this section only discusses one-dimensional genetic
values.

```{eval-rst}
.. todo:: 

    Update for multi-dimensional genetic values once we get documentation
    in place for simulations of Gaussian stabilizing selection with pleiotropy.
```

(python_genetic_value_update_desc)=

#### Supporting time-dependent models

Genetic value methods may have dynamic states that depend on some
state of the population.  For example, something may change once
`pop.generation > x`.

To support such time-dependent states, classes have a method with
the following signature:

```{code-block} python

def update(self, pop: fwdpy11.DiploidPopulation) -> None:
    pass

```

A type with no time-dependent state can in fact use the function shown
above. A concrete example is shown {ref}`below <python-gvalue-to-fitness>`,
in an example class called `PyGSSRandomOptimum`.  That example does not
use any information from `pop` to update its state.

### Performance

It will be hard to get "awesome" performance with Python genetic
value classes.  The culprit are abstract methods that get called
once per offspring.  That's simply a lot of Python/C++ boundary
crossing.  The `update` function described above is called only
once per generation, making it much less of a performance killer.

We've tried hard to make the data exchange from C++ to your Python
code as efficient as possible.  We accomplished this by carefully
considering how we implement the `data` argument that gets passed
in.  Benchmarking showed that carefully designed proxies to the
underlying C++ data greatly improved run times.

Importantly, you do not have to write every component in Python!
If you are interested in additive effects models but a new "noise"
model, then use {class}`fwdpy11.Additive` along with your custom
type.

`fwdpy11` provides the following helper functions to improve efficiency:

```{py:function} fwdpy11.strict_additive_effects

This function computes the sum of effect sizes in a
diploid genome, ignoring dominance.  The return value
is an offset from zero, so add 1.0 to convert to a
fitness if needed.

:param pop: The population
:type pop: fwdpy11.DiploidPopulation
:param metadata: The offspring metadata
:type metadata: fwdpy11.DiploidMetadata
:returns: Sum of effect sizes
:rtype: float

```

```{py:function} fwdpy11.additive_effects

This function computes the sum of effect sizes in a
diploid genome, accounting for dominance.  The return value
is an offset from zero, so add 1.0 to convert to a
fitness if needed.

:param pop: The population
:type pop: fwdpy11.DiploidPopulation
:param metadata: The offspring metadata
:type metadata: fwdpy11.DiploidMetadata
:type scaling: float
:returns: Sum of effect sizes (accounting for heterozygous effects)
:rtype: float

```

The first function is equivalent to the following Python code:

```{code-block} python

def strict_additive_effects(
    pop: fwdpy11.DiploidPopulation,
    metadata: fwdpy11.DiploidMetadata,
) -> float:
    g = 0.0
    dip = pop.diploids[metadata.label]
    for i in [dip.first, dip.second]:
        for k in pop.haploid_genomes[i].smutations:
            g += pop.mutations[k].s
    return g

```

The second function uses C++ code from `fwdpp` to do the calculation
accounting for dominance effects.

### Technical issues

This section discusses some important technicalities.

#### Abstract base classes

In object oriented lingo, the `fwdpy11` base classes here are *abstract
base classes* or `ABCs`.  This means that they simply describe a required
interface.  When you inherit from an ABC, you must fill in the details
of that interface.

#### Initializing the super class

When writing object-oriented code, it is important that a derived
class properly initialize the base class.  In Python, base
classes are often called `superclasses` and derived classes are
`subclasses`.

The idiomatic approach in `Python 3` is:

```{code-cell} python
class Base(object):
   def __init__(self):
       self.b = "I am a Base"
```

```{code-cell} python
class Derived(Base):
   def __init__(self):
       self.d = "I am a Derived"
       super(Derived, self).__init__()
```

```{code-cell} python
d = Derived()
print(d.b)
```

Failure to do that `super(...)` business would result in an
`AttributeError` when trying to print `d.b`.

However, the classes that you are writing are **not** straightforward
Python classes!  While `super(...)` *may* work, it is not
guaranteed.  To quote the `pybind11` documentation:

> Note that a direct __init__ constructor should be called,
> and super() should not be used. For simple cases of linear inheritance,
> super() may work, but once you begin mixing Python and C++ multiple inheritance,
> things will fall apart due to differences between Python’s MRO and C++’s mechanisms.
> 
The source of this quote is the `pybind11` manual section on
[Classes](<https://pybind11.readthedocs.io/en/stable/advanced/classes.html#>).

The examples that follow take their advice.

#### `attrs`

The [attrs](<https://www.attrs.org>) package provides a library
of convenient class decorators that significantly reduce the amount
of "boiler plate" code needed to write Python
classes.  Some of the examples below will use this library.  You must
*not* use `slots=True`!  The `superclasses` already use the slots mechanism,
and it seems difficult for a derived class to also use slots and properly
initialize the base.

Other than this caveat regarding slots, we enthusiastically recommend
`attrs`.

#### Custom decorators

Speaking of decorators, we provide several that make writing genetic
value types a bit easier.  You will see concrete examples using them
later.

```{py:decorator} fwdpy11.custom_genetic_value_decorators.genetic_value_is_trait_default_clone(ndim=1)

This decorator is a *class* that adds the ability
for custom types to be passed to classes derived from
:class:`fwdpy11.DiploidGeneticValue`.

Because this decorator is a class, you need the ``()``
to apply it.

This decorator is required because of the very different
object models of C++ and Python.
See `here <https://github.com/pybind/pybind11/issues/1049>`_
for some of the gory details.

.. py:method:: __init__

    :param ndim: The dimensionality of the model
    :type ndim: int

    The default is ``ndim=1``.

```

```{py:decorator} fwdpy11.custom_genetic_value_decorators.genetic_value_noise_default_clone

This is the analog to
:attr:`fwdpy11.custom_genetic_value_decorators.genetic_value_is_trait_default_clone`
for noise objects inheriting from :class:`fwdpy11.GeneticValueNoise`.

Unlike its analog, this decorator is a *function* and not a class.

```

```{py:decorator} fwdpy11.custom_genetic_value_decorators.default_update

This decorator adds the following function to ``cls``,
which is a no-op function:

.. code-block:: python

    def update(pop: fwdpy11.DiploidPopulation):
        pass

.. note::

    In the future, we will experiment with adding a C++
    implementation of this no-op.

```

(python-gvalue-to-fitness)=

### Genetic value to fitness maps

Here, one derives a new class from {class}`fwdpy11.GeneticValueIsTrait`.
The basic form of such a class is:

```{code-block} python

import fwdpy11
import fwdpy11.custom_genetic_value_decorators


@fwdpy11.custom_genetic_value_decorators.genetic_value_is_trait_default_clone()
class MyGeneticValueToFitness(fwdpy11.GeneticValueIsTrait):
    def __init__(self):
        # Do NOT call super()!
        fwdpy11.GeneticValueIsTrait.__init__(self)

    def __call__(self, data: fwdpy11.DiploidGeneticValueToFitnessData) -> float:
        pass

    def update(pop: fwdpy11.DiploidPopulation) -> None:
        pass

```

If the object does not need to manage an internal state that changes over time,
we can simplify a bit with another decorator:

```{code-block} python

import fwdpy11
import fwdpy11.custom_genetic_value_decorators


@fwdpy11.custom_genetic_value_decorators.default_update
@fwdpy11.custom_genetic_value_decorators.genetic_value_is_trait_default_clone()
class MyGeneticValueToFitness(fwdpy11.GeneticValueIsTrait):
    def __init__(self):
        # Do NOT call super()!
        fwdpy11.GeneticValueIsTrait.__init__(self)

    def __call__(self, data: fwdpy11.DiploidGeneticValueToFitnessData) -> float:
        pass

```

Two complete examples come from the test suite.  The first example reimplements
{class}`fwdpy11.GSS` in Python.  The second example implements a model
where the optimum changes to a new value each generation. (Note that
the second model requires that {func}`numpy.random.seed` be called elsewhere
so that results are reproducible.)

```{literalinclude} ../../tests/pygss.py
:lines: 19-

```

#### The data type

The type passed into the `__call__` function is:

```{py:class} fwdpy11.DiploidGeneticValueToFitnessData

.. versionadded:: 0.9.0

This class supports the buffer protocol, which exposes the
genetic values array.  The most efficient access will
be via a :class:`memoryview`.

Instances of this class have the following attributes:

.. py:attribute:: offspring_metadata

    An instance of :class:`fwdpy11.DiploidMetadata`, giving you
    access to those fields that have been assigned to the offspring.
    (This is a *copy* of the metadata from the C++ side.)

.. py:attribute:: offspring_metadata_index

    A 64 bit integer that gives the location (index) of
    ``offspring_metadata`` in :attr:`fwdpy11.DiploidPopulation.diploid_metadata`.
    This index is useful in the event of mass migrations via copies,
    which can cause a mismatch between :attr:`fwdpy11.DiploidMetadata.label`
    and this value.

.. py:attribute:: parent1_metadata

    The first parent's metadata (:class:`fwdpy11.DiploidMetadata`).

.. py:attribute:: parent2_metadata

    The second parent's metadata (:class:`fwdpy11.DiploidMetadata`).

```

### Genetic value "noise"

A custom "noise" class inherits from {class}`fwdpy11.GeneticValueNoise`.
A minimal implementation has the following form:

```{code-block} python

@fwdpy11.custom_genetic_value_decorators.genetic_value_noise_default_clone
class MyGeneticValueNoise(fwdpy11.GeneticValueNoise):
    def __call__(self, data: fwdpy11.DiploidGeneticValueNoiseData) -> float:
        pass

    def update(pop: fwdpy11.DiploidPopulation) -> None:
        pass

```

If your model can allow a no-op `update` function, you may apply
the decorator
{attr}`fwdpy11.custom_genetic_value_decorators.genetic_value_noise_default_clone`.

A simple example from the test suite implements Gaussian noise with mean zero
and standard deviation `0.1`:

```{literalinclude} ../../tests/pynoise.py
:lines: 19-

```

See {func}`fwdpy11.gsl_ran_gaussian_ziggurat` for details on that function.

#### The data type

```{py:class} fwdpy11.DiploidGeneticValueNoiseData

.. versionadded:: 0.9.0

Instances of this class have the following attributes:

.. py:attribute:: rng

    The simulation's random number generation, an
    instance of :class:`fwdpy11.GSLrng`

.. py:attribute:: offspring_metadata

    An instance of :class:`fwdpy11.DiploidMetadata`, giving you
    access to those fields that have been assigned to the offspring.
    (This is a *copy* of the metadata from the C++ side.)

.. py:attribute:: offspring_metadata_index

    A 64 bit integer that gives the location (index) of
    ``offspring_metadata`` in :attr:`fwdpy11.DiploidPopulation.diploid_metadata`.
    This index is useful in the event of mass migrations via copies,
    which can cause a mismatch between :attr:`fwdpy11.DiploidMetadata.label`
    and this value.

.. py:attribute:: parent1_metadata

    The first parent's metadata (:class:`fwdpy11.DiploidMetadata`).

.. py:attribute:: parent2_metadata

    The first parent's metadata (:class:`fwdpy11.DiploidMetadata`).

```

### New genetic value models

#### The abstract base class

The ABC type is {class}`fwdpy11.PyDiploidGeneticValue`:

```{py:class} fwdpy11.PyDiploidGeneticValue

.. py:method:: __init__()

    :param ndim:
    :type ndim: int
    :param genetic_value_to_fitness:
    :type genetic_value_to_fitness: fwdpy11.GeneticValueIsTrait or None
    :param noise:
    :type noise: fwdpy11.GeneticValueNoise or None

.. py:method:: calculate_gvalue
    :abstractmethod:

    :param data: Input data
    :type data: fwdpy11.PyDiploidGeneticValueData
    :returns: The value to be stored in the offspring's
              :attr:`fwdpy11.DiploidMetadata.g`
    :rtype: float

.. py:method:: genetic_value_to_fitness

    :param data: Input data
    :type data: fwdpy11.DiploidGeneticValueToFitnessData

    :returns: fitness
    :rtype: float

    .. note::

        This function does not need to be defined
        by derived classes most of the time.
        The default behavior is to apply
        the instance of :class:`fwdpy11.GeneticValueToFitnessMap`
        stored by the instance of :class:`fwdpy11.DiploidGeneticValue`
        (or :class:`fwdpy11.PyDiploidGeneticValue`).  Defining
        this function in a derived class skips calling held instance
        in favor of the derived class implementation.
        In general, one only needs to derive this class for models
        where either individual genetic values depend on genotypes of the
        rest of the population. See :ref:`here <more-complex-gvalue-models>`.


.. py:method:: update

    :param pop: The population being simulated
    :type pop: fwdpy11.DiploidPopulation
    :rtype: None

    .. note::

        A default implementation can be defined using
        the decorator
        :attr:`fwdpy11.custom_genetic_value_decorators.genetic_value_noise_default_clone`.

```

#### Outline of user-defined classes

The form of a custom genetic value type inherits from the ABC
described above and would look something like this:

```{code-block} python

class MyGeneticValueObject(fwdpy11.PyDiploidGeneticValue):
    def __init__(self, *args):
        # Do not call super()!
        fwdpy11.PyDiploidGeneticValue.__init__(...)

    def calculate_gvalue(self, data: fwdpy11.PyDiploidGeneticValueData) -> float:
        pass

    def update(self, pop: fwdpy11.DiploidPopulation) -> None:
        pass

```

#### Method requirements

{class}`fwdpy11.PyDiploidGeneticValue` *does* provide a default `update`
function that updates the genetic-value-to-fitness-map and noise objects.

If your model can allow a no-op `update` function, you may apply
the decorator
{attr}`fwdpy11.custom_genetic_value_decorators.genetic_value_noise_default_clone`.

The `genetic_value_to_fitness` function must completely replace the
functionality of an instance of {class}`fwdpy11.GeneticValueIsTrait`.
See {ref}`here <python-gvalue-to-fitness>` for details.

The `calculate_gvalue` function is more complex because it has more
responsibilities:

1. It must return a value that will populate the genetic value
   field of the offspring's metadata ({attr}`fwdpy11.DiploidMetadata.g`).
2. It also needs to populate an array of {class}`float` that represents
   the genetic values for all trait dimensions. See {ref}`here <python_genetic_value_dimensions>`.
   For the common case of a one-dimensional genetic value,
   you populate element `0` of this array with value that you
   will return from this function. The "genetic values array" is
   accessible via the Python buffer protocol.
   The next section shows a concrete example.

:::{note}

The steps described in `2` are populating
{attr}`fwdpy11.DiploidGeneticValue.genetic_values`.  In the
future, we may provide read-write access to that field, which
would constitute an API change.  However, the path of least
resistance in the short term was to implement the current approach.

:::

#### Example: strictly additive effects on a trait

Here is an example of a strict additive model of a trait,
taken from the test suite. You'll notice that we are accessing
arrays on the `C++` using Python {class}`memoryview` objects.
This is a *very* efficient method, and much faster
than creating a temporary, non-owning {class}`numpy.ndarray`.

```{literalinclude} ../../tests/pyadditive.py
:lines: 20-

```

If you want to account for the heterozygous effect of mutations, then the
most explicit approach is probably to use {class}`collections.Counter`:

```{code-block} python

def calculate_gvalue(self, data):
    s = 0.0
    c = Counter()
    pop = data.pop
    dip = data.diploids[data.offspring_metadata.label]
    for g in [dip.first, dip.second]:
        for k in pop.haploid_genomes[g].smutations:
            m = pop.mutations[k]
            c.update([(m.pos, m.s, m.h)])
    for i, j in c.items():
        if j == 1:  # Aa
            s += i[1] * i[2]
        else:  # aa
            # This is a hard-coded "scaling"
            # of 2
            s += 2.0 * i[1]
    # Use a memory view to update
    # the genetic values array
    memoryview(data)[0] = s
    return s

```

For additive cases like this, though, you'll get better performance
with {func}`fwdpy11.additive_effects`.

(more-complex-gvalue-models)=

#### More complex scenarios (such as social interactions)

Some models are not easily implemented with the genetic value +
genetic value to fitness map + noise trio of types defined here.
The reason has to do with how these types are applied during a simulation,
which is described in more detail {ref}`here <gvalue-technical-details>`, and
in {ref}`this section <gvalue-simulation-behavior>` specifically.

#### The data types

```{eval-rst}
.. class:: fwdpy11.PyDiploidGeneticValueData

.. versionadded:: 0.9.0

This class supports the buffer protocol, which exposes the
genetic values array.  The most efficient access will
be via a :class:`memoryview`.

Instances of this class have the following attributes:

.. py:attribute:: rng

    The simulation's random number generation, an
    instance of :class:`fwdpy11.GSLrng`

.. py:attribute:: pop

    The population, an instance of :class:`fwdpy11.DiploidPopulation`

.. py:attribute:: offspring_metadata

    The offspring's :class:`fwdpy11.DiploidMetadata`

.. py:attribute:: parent1_metadata

    The first parent's :class:`fwdpy11.DiploidMetadata`

.. py:attribute:: parent2_metadata

    The second parent's :class:`fwdpy11.DiploidMetadata`

.. py:attribute:: offspring_metadata_index

    A 64 bit integer that gives the location (index) of
    ``offspring_metadata`` in :attr:`fwdpy11.DiploidPopulation.diploid_metadata`.
    This index is useful in the event of mass migrations via copies,
    which can cause a mismatch between :attr:`fwdpy11.DiploidMetadata.label`
    and this value.
```

(msprime-subtleties)=

## Subtleties of starting/finishing tree sequences using `msprime`

Some care is required starting a simulation from a tree sequence (see {ref}`here <precapitation>`) or finishing a simulation with one (see {ref}`here <recapitation>`).

If simulating with discrete recombination positions (see {ref}`here <geneticmapunit>`), it may be important to realize that the current stable releases of `msprime` model continuous recombination breakpoints.
A future `msprime` release will support discretized positions.

When simulating with genome lengths longer than unity (1.0), determining the correct recombination rate for `msprime` needs additional work.
The examples shown {ref}`here <recapitation>` and {ref}`here <precapitation>` show the following method for determining the `msprime` recombination rate:

```{code-block} python

recombination_rate = rho / Ne / 4

```

That code implies a genome length of 1.0.
A more general example, for a genome of length `L` would look like:

```{code-block} python

recombination_rate = rho / Ne / 4 / L

```


