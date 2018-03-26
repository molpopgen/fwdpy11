.. _genetic_values_types:

Objects for calculation of genetic values
====================================================================================

Background reading:

* :ref:`genetic_values`

For the calculation of *fitness* in a single-region simulation, the following objects are available:

* :class:`fwdpy11.fitness.SlocusAdditive`
* :class:`fwdpy11.fitness.SlocusMult`

These two classes take a **scaling** parameter as a constructor argument:

.. ipython:: python

    import fwdpy11.fitness

    # For a single mutation, fitness is
    # 1, 1+sh, 1+s
    a1 = fwdpy11.fitness.SlocusAdditive(1.0)

    # For a single mutation, fitness is
    # 1, 1+sh, 1+2*s
    a2 = fwdpy11.fitness.SlocusAdditive(2.0)

For the calculation of *genetic values* in a single-region simulation, the following objects are available:

* :class:`fwdpy11.trait_values.SlocusAdditiveTrait`
* :class:`fwdpy11.trait_values.SlocusMultTrait`
* :class:`fwdpy11.trait_values.SlocusGBRTrait`

The first two objects also take **scaling** parameters as constructor arguments.  The third class is only valid for
distributions of effects sizes returning non-negative values.

For multi-region simulations, calculations are run through instances of
:class:`fwdpy11.multilocus.MultiLocusGeneticValue`.  Instances of this object hold lists of single-region objects:

.. ipython:: python

    import fwdpy11.trait_values
    import fwdpy11.multilocus

    # Genetic value calculator for a quantitative trait simulation
    # of ten regions under a "Turelli-like" additive model of 
    # genotype -> phenotype:
    am = fwdpy11.multilocus.MultiLocusGeneticValue([fwdpy11.trait_values.SlocusAdditiveTrait(2.0)]*10)

During a simulation, each function stored is applied to each region.  The result is a list of floating-point values,
which are the genetic values due to each locus.  These have to be *aggregated* into a final genetic value for an
individual.  An aggregator is any callable accepting a NumPy array and returning a single float.  We provide the
following built-in aggregator types:

* :class:`fwdpy11.multilocus.AggAddFitness`
* :class:`fwdpy11.multilocus.AggMultFitness`
* :class:`fwdpy11.multilocus.AggAddTrait`
* :class:`fwdpy11.multilocus.AggMultTrait`

The first two are for *fitness* calculations and the latter two for *genetic value* calculations (for quantitative trait
simulations).

You may provide your own aggregators, written either in C++ or in Python.

.. note::

    It is really important to match your aggregator types to your single-locus types.  Mixing a list of
    :class:`fwdpy11.trait_values.SlocusAdditiveTrait` with :class:`fwdpy11.multilocus.AggAddFitness` for
    an aggregator will give strange results, as you are mixing a zero-centered calculator with an aggregator centered on
    one.

The relationship to fixations
--------------------------------------------------------------------

For standard population-genetic simulations, relative fitness is what matters.  Relative fitnesses are unaffected by
fixations under multiplicative models, but the same is not true under additive models.  Please note that multiplicative
models are typically assumed, and thus you should use :class:`fwdpy11.fitness.SlocusMult` most of the time.  Doing so
will simply make your life easier (and your simulations more efficient--keep reading...).

For simulations of phenotypes where fitness is determined by comparing phenotype to some optimum value, fixations always
affect the distance of an individual from this optimum.

The reason to bring all this up is because fixations may be removed from gametes during simulation, depending on
parameters that you input.  Pruning fixations results in faster simulations, because those sites are not considered in
fitness calculations.  However, you should *not* prune them when simulating additive models of fitness or when
simulating phenotypes.  See :ref:`handling_fixations` for more details.

The future
-----------------------------------------------------------

We hope to:

* Reduce the number of single locus types so that the trait-vs-fitness decision is a constructor argument.
* Add a GBR type for fitness.
* Reduce the complexity of the code underyling all this, which will hopefully enable the above.
* Make all this stuff about fixations something that the user (you) doesn't have to worry about.
