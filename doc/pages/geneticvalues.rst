.. _genetic_values_types:

Objects for calculation of genetic values
====================================================================================

.. versionchanged:: 0.1.5

    Updated to cover new types in :py:mod:`fwdpy11.genetic_values`. These type unify several concepts into 
    a simpler Python framework with a simpler C++ back-end.

Background reading:

* :ref:`genetic_values`

In order to calculate a diploid's fitness, we need to process its mutations, account for any random effects, and produce
a final fitness value.  In other words, a diploid's _genetic value_ is a function of its mutations, and the final
fitness is the result of mapping that genetic value to fitness (accounting for any random effects). 

Depending on the type of simulation we are doing, a diploid's genetic value, :math:`G`, may represent fitness
(:math:`w`) itself or it may represent a trait value (a "phenotype").  The module :py:mod:`fwdpy11.genetic_values` provides a flexible Python
class hierarchy to accomodate these different scenarios.  There are two different class hiarchies.  One is used for
simulations involving :class:`fwdpy11.SlocusPop` and the other applies to :class:`fwdpy11.MlocusPop`.  The base class of
each is an Abstract Base Class (ABC):

* :class:`fwdpy11.genetic_values.SlocusPopGeneticValue`
* :class:`fwdpy11.genetic_values.MlocusPopGeneticValue`

Because these two classes are ABCs, you may not make instances of them.  They exist to provide the following minimal
interface to the user:

* They are callable types, capable of returning the genetic value of the :math:`i^{th}` diploid in a :class:`fwdpy11.Population`.
  See the documentation for :func:`fwdpy11.genetic_values.SlocusPopGeneticValue.__call__` and
  :func:`fwdpy11.genetic_values.MlocusPopGeneticValue.__call__`.
* The may return the fitness of the :math:`i^{th}` diploid.  See docstrings for :func:`fwdpy11.genetic_values.SlocusPopGeneticValue.fitness`
  and :func:`fwdpy11.genetic_values.MlocusPopGeneticValue.fitness`.

Internally, objects in these class hierarchies provide the following functionality:

* Calculate genetic value based on mutations.
* Calculate any random effects, or "noise".
* Map genetic value to fitness.

Two more ABCs define Python classes capable of flexibly modeling noise and mappings of genetic values to fitness:

* :class:`fwdpy11.genetic_values.SlocusPopGeneticValueWithMapping`
* :class:`fwdpy11.genetic_values.MlocusPopGeneticValueWithMapping`

These two types inherit from the two ABCs described above, and thus provide the same public interface.  They
additionally provide:

* The ability to access the genetic value to fitness map.  
* Acccess to the "noise" function.

At this point, it helps to look at a concrete class, which will allow us to look at the behavior of the interfaces
defined by the ABCs:

.. ipython:: python

    import fwdpy11.genetic_values

    multiplicative = fwdpy11.genetic_values.SlocusMult(2.0)

    print(type(multiplicative.gvalue_to_fitness))
    print(type(multiplicative.noise))

In the above code, we created an instance of :class:`fwdpy11.genetic_values.SlocusMult`, which models multiplicative
genetic values.  Our mapping of genetic value to fitness is handled by an instance of
:class:`fwdpy11.genetic_values.GeneticValueIsFitness`.  As the name implies, we will be simulating mutations with
_direct_ effects on fitness.  Thus, in the absence of random effects, :math:`w = G`.  Here, the type generating random
effets on genetic values is :class:`fwdpy11.genetic_value_noise.NoNoise`.  Again, the name should make it obvious what is
going on: there are no random effects!  Thus, the variable `multiplicative` will model the standard population genetic
scenario of multiplicative mutational effects on fitness.

.. note::

    You have just learned that the types handling noise are in :py:mod:`fwdpy11.genetic_value_noise`.

Let's look at a few more properties of our variable:

.. ipython:: python

    # Does this type model fitness or a trait?
    print(multiplicative.is_fitness)
    # What is the scaling parameter?
    print(multiplicative.scaling)
    # What are the relations of this type to our class hierarchy?
    print(isinstance(multiplicative, fwdpy11.genetic_values.SlocusPopGeneticValueWithMapping))
    print(isinstance(multiplicative, fwdpy11.genetic_values.SlocusPopGeneticValue))

.. note::

    For more details on the scaling parameter, see the documentation for :class:`fwdpy11.genetic_values.SlocusMult`.

Let's look at an example where :math:`G \neq w` and there are random effects:

.. ipython:: python

    import fwdpy11.genetic_value_noise

    mult_trait = fwdpy11.genetic_values.SlocusMult(2.0, 
        fwdpy11.genetic_values.GSS(opt = 0.0, VS = 1.0),
        fwdpy11.genetic_value_noise.GaussianNoise(mean=0.0, sd=0.1))
    print(mult_trait.is_fitness)

Now, we have a model where :math:`G` is a genetic value determined by multiplicative interactions amongst mutations.
Gaussian noise with mean zero and standard deviation 0.1 is added to :math:`G` to determine the final phenotype and
fitness, :math:`w` is modeled by Gaussian stabilizing selection (GSS) with an optimum trait value of zero and a strength
of stabilizing selection, :math:`VS`, equal to one.

The following types are provided in :py:mod:`fwdpy11.genetic_values` to calculate genetic_values/fitness:

* :class:`fwdpy11.genetic_values.SlocusMult`
* :class:`fwdpy11.genetic_values.SlocusAdditive`
* :class:`fwdpy11.genetic_values.SlocusGBR`
* :class:`fwdpy11.genetic_values.MlocusMult`
* :class:`fwdpy11.genetic_values.MlocusAdditive`
* :class:`fwdpy11.genetic_values.MlocusGBR`

.. note::

    The "GBR" types are only usable for models of quantitative traits and not for models of direct effects on fitness.
    See the relevant papers, which are cited in the docstrings for the classes, for details.

The following types map genetic value to fitness:

* :class:`fwdpy11.genetic_values.GeneticValueIsFitness`
* :class:`fwdpy11.genetic_values.GSS`
* :class:`fwdpy11.genetic_values.GSSmo`

In the above list, the first type (:class:`fwdpy11.genetic_values.GeneticValueIsFitness`) is used for "standard
population genetic" simulations.  In our first example code block above, we see that it is used as a default value.  The
latter two classes are used to model quantitative traits.

To learn more about random effects, see :py:mod:`fwdpy11.genetic_value_noise`.


The relationship to fixations
--------------------------------------------------------------------

For standard population-genetic simulations, relative fitness is what matters.  Relative fitnesses are unaffected by
fixations under multiplicative models, but the same is not true under additive models.  Please note that multiplicative
models are typically assumed, and thus you should use :class:`fwdpy11.genetic_values.SlocusMult` or :class:`fwdpy11.genetic_values.MlocusMult` 
most of the time.  Doing so will simply make your life easier (and your simulations more efficient--keep reading...).

For simulations of phenotypes where fitness is determined by comparing phenotype to some optimum value, fixations always
affect the distance of an individual from this optimum.

The reason to bring all this up is because fixations may be removed from gametes during simulation, depending on
parameters that you input.  Pruning fixations results in faster simulations, because those sites are not considered in
fitness calculations.  However, you should *not* prune them when simulating additive models of fitness or when
simulating phenotypes.  See :ref:`handling_fixations` for more details.

The future
-----------------------------------------------------------

We hope to:

* Add a GBR type for fitness.
* Make all this stuff about fixations something that the user (you) doesn't have to worry about.
