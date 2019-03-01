.. _popobjects:

Construction of population objects from a predetermined state
============================================================================================================================================

Background reading:

* :ref:`data_types`

Beginning with version 0.1.4, you are able to construct population objects from pre-calculated containers of diploid,
gamete, and mutation objects.  Constructing population objects this way allows for some very powerful simulation
methods.  However, some care is required, and you should be familiar with the material in :ref:`data_types`, which
covers how populations are represented in memory.

.. note::
    This page uses :class:`fwdpy11.SlocusPop` in its examples.
    However, analagous functionality is available for all supported population types.
    
Quick start example
-----------------------------------

Building populations this way will take advantage of object constructors taking tuples as parameters.  Let's just dive
in:

.. testcode:: constructing_pops

    import fwdpy11

    # Create empty containers for mutations,
    # gametes, and diploids:
    mutations = fwdpy11.VecMutation()
    gametes = fwdpy11.VecGamete()
    diploids = fwdpy11.VecDiploid()

    # Add a mutation with pos = 0.1,
    # s = -0.01, h = 1.0, g = 0,
    # and label = 0
    mutations.append(fwdpy11.Mutation(0.1,-0.01,1.0,0,0))

    # Add a gamete that exists in two copies
    # and contains our mutation.  The mutation
    # is not neutral, and therefore belongs in
    # 'smutations', with is the second vector
    # in our tuple.  The first is the list of 
    # indexes of neutral variants.
    gametes.append(fwdpy11.Gamete((2,fwdpy11.VecUint32([]),fwdpy11.VecUint32([0]))))

    # There is only one gamete, and so our diploid
    # will contain two copies of it, and its index is 
    # 0:
    diploids.append(fwdpy11.DiploidGenotype(0,0))
    
    # Create a population from our containers
    pop = fwdpy11.SlocusPop(diploids, gametes, mutations)

    # Check that our pop is as expected:
    assert(pop.N == 1)
    assert(len(pop.diploids) == 1)
    assert(len(pop.gametes) == 1)
    assert(len(pop.mutations) == 1)
    assert(pop.gametes[0].n == 2)
    assert(mutations[0].pos == 0.1)
    assert(mutations[0].s == -0.01)
    assert(mutations[0].h == 1)
    assert(mutations[0].g == 0)
    assert(mutations[0].label == 0)


The above example brings up some new concepts:

1. There are several new container types introduced. These are thin wrappers around C++ containers. 
2. We _can_ construct these containers from Python iterable types such as lists.  The above example does just that to
   create a :class:`fwdpy11.Gamete` instance.  However, doing so with large lists may use a lot of RAM, and
   the append approach may be preferable.

Some comments are needed:

1. The C++ back-end (fwdpp_) is very strict about the input.  Any error will result in exceptions being thrown.

Efficient construction of large populations
-----------------------------------------------

When building a large population programmatically, a naive approach would leave us with a data copy on both the Python
and on the C++ side, which is not ideal.  We can avoid that via class methods that use C++11's move semantics to "steal"
the data from our input containers. The following example is the same as above, except that we create the function via
:func:`fwdpy11.SlocusPop.create`:

.. testcode:: move_constructing_pops

    import fwdpy11

    mutations = fwdpy11.VecMutation()
    gametes = fwdpy11.VecGamete()
    diploids = fwdpy11.VecDiploid()

    mutations.append(fwdpy11.Mutation(0.1,-0.01,1.0,0,0))

    gametes.append(fwdpy11.Gamete((2,fwdpy11.VecUint32([]),fwdpy11.VecUint32([0]))))

    diploids.append(fwdpy11.DiploidGenotype(0,0))
    
    pop = fwdpy11.SlocusPop.create(diploids, gametes, mutations)
    assert(len(diploids) == 0)
    assert(len(gametes) == 0)
    assert(len(mutations) == 0)
    assert(len(pop.diploids) == 1)
    assert(len(pop.mutations) == 1)
    assert(len(pop.gametes) == 1)

The first three assertions show that the containers that we contstructed are now empty.  Their contents have been moved
into the population object, avoiding an extra temporary copy.  *The create function should be the preferred method of
constructing populations unless you have a reason to keep the input data around.*


Examples of input errors
-----------------------------------------------

Incorrect gamete count:

.. testcode::

    import fwdpy11
    mutations = fwdpy11.VecMutation()
    gametes = fwdpy11.VecGamete()
    diploids = fwdpy11.VecDiploid()
    mutations.append(fwdpy11.Mutation(0.1,-0.01,1.0,0,0))
    # The gamete is incorrectly labelled as occurring once:
    gametes.append(fwdpy11.Gamete((1,fwdpy11.VecUint32([]),fwdpy11.VecUint32([0]))))
    diploids.append(fwdpy11.DiploidGenotype(0,0))
    pop = fwdpy11.SlocusPop.create(diploids, gametes, mutations)

The result is a `RuntimeError`:

.. testoutput::
    :options: +ELLIPSIS

    Traceback (most recent call last):
    ...
    RuntimeError: gamete count does not match number of diploids referring to it

Neutral or non-neutral mutations in the incorrect gamete container:

.. testcode::

    import fwdpy11
    mutations = fwdpy11.VecMutation()
    gametes = fwdpy11.VecGamete()
    diploids = fwdpy11.VecDiploid()
    mutations.append(fwdpy11.Mutation(0.1,-0.01,1.0,0,0))
    # The mutation is non-neutral, and we are mistakenly
    # putting it in the Gametes.mutations container:
    gametes.append(fwdpy11.Gamete((2,fwdpy11.VecUint32([0]),fwdpy11.VecUint32([]))))
    diploids.append(fwdpy11.DiploidGenotype(0,0))
    pop = fwdpy11.SlocusPop.create(diploids, gametes, mutations)

The result is a `RuntimeError`:

.. testoutput::
    :options: +ELLIPSIS

    Traceback (most recent call last):
    ...
    RuntimeError: gamete contains key to mutation in wrong container.

Other conditions that will lead to errors include:

1. Gametes and diploids containing indexes that are out of range.
2. Mutation keys in gametes must be sorted according to mutation position.


Dealing with unsorted mutation input 
---------------------------------------------------------------------------------------------------------
Consider the following example with two mutations:

.. ipython:: python

    import fwdpy11

    mutations = fwdpy11.VecMutation()

    gametes = fwdpy11.VecGamete()

    diploids = fwdpy11.VecDiploid()

    mutations.append(fwdpy11.Mutation(0.1,-0.01,1.0,0,0))

Add in a second, non-neutral mutation:

.. ipython:: python

    mutations.append(fwdpy11.Mutation(0.22,0.1,1.0,0,1))

Put mutations into containers out of order as far as mutation position is concerned:

.. ipython:: python

    gametes.append(fwdpy11.Gamete((2,fwdpy11.VecUint32([]),fwdpy11.VecUint32([1,0]))))

    diploids.append(fwdpy11.DiploidGenotype(0,0))

We will get an exception when we try to create a population:

.. ipython:: python
    :okexcept:

    pop = fwdpy11.SlocusPop(diploids, gametes, mutations)

We can sort the input data with a call to :func:`fwdpy11.util.sort_gamete_keys`:

.. ipython:: python

    from fwdpy11.util import sort_gamete_keys
    sort_gamete_keys(gametes,mutations)
    pop = fwdpy11.SlocusPop.create(diploids, gametes, mutations)

The sorting takes place on the C++ side because of how the relevant container types are exposed to Python.

.. testcode::
    :hide:

    import fwdpy11
    mutations = fwdpy11.VecMutation()
    gametes = fwdpy11.VecGamete()
    diploids = fwdpy11.VecDiploid()
    mutations.append(fwdpy11.Mutation(0.1,-0.01,1.0,0,0))
    # Add in a second, non-neutral mutation:
    mutations.append(fwdpy11.Mutation(0.22,0.1,1.0,0,1))
    # Put mutations into containers out of order
    # as far as mutation position is concerned:
    gametes.append(fwdpy11.Gamete((2,fwdpy11.VecUint32([]),fwdpy11.VecUint32([1,0]))))
    diploids.append(fwdpy11.DiploidGenotype(0,0))
    pop = fwdpy11.SlocusPop.create(diploids, gametes, mutations)

.. testoutput::
    :hide:

    Traceback (most recent call last):
    ...
    ValueError: gamete contains unsorted keys 

.. _fwdpp: http://molpopgen.github.io/fwdpp
.. _work: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3969567/
