.. _popobjects:

Construction of population objects from a predetermined state.
============================================================================================================================================

Background reading:

* :ref:`data_types`

Beginning with version 0.1.4, you are able to construct population objects from pre-calculated containers of diploid,
gamete, and mutation objects.  Constructing population objects this way allows for some very powerful simulation
methods.  However, some care is required, and you should be familiar with the material in :ref:`data_types`, which
covers how populations are represented in memory.

Quick start example
-----------------------------------

Building populations this way will take advantage of object constructors taking tuples as parameters.  Let's just dive
in:

.. testcode:: constructing_pops

    import fwdpy11

    # Create empty containers for mutations,
    # gametes, and diploids:
    mutations = fwdpy11.MutationContainer()
    gametes = fwdpy11.GameteContainer()
    diploids = fwdpy11.DiploidContainer()

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
    gametes.append(fwdpy11.Gamete((2,fwdpy11.VectorUint32([]),fwdpy11.VectorUint32([0]))))

    # There is only one gamete, and so our diploid
    # will contain two copies of it, and its index is 
    # 0:
    diploids.append(fwdpy11.SingleLocusDiploid(0,0))
    
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
   create a :class:`fwdpy11.fwdpp_types.Gamete` instance.  However, doing so with large lists may use a lot of RAM, and
   the append approach may be preferable.

Some comments are needed:

1. The C++ back-end (fwdpp_) is very strict about the input.  Any error will result in exceptions being thrown.

Examples of input errors
-----------------------------------------------

Incorrect gamete count:

.. testcode::

    import fwdpy11
    mutations = fwdpy11.MutationContainer()
    gametes = fwdpy11.GameteContainer()
    diploids = fwdpy11.DiploidContainer()
    mutations.append(fwdpy11.Mutation(0.1,-0.01,1.0,0,0))
    # The gamete is incorrectly labelled as occurring once:
    gametes.append(fwdpy11.Gamete((1,fwdpy11.VectorUint32([]),fwdpy11.VectorUint32([0]))))
    diploids.append(fwdpy11.SingleLocusDiploid(0,0))
    pop = fwdpy11.SlocusPop(diploids, gametes, mutations)

The result is a `RuntimeError`:

.. testoutput::
    :options: +ELLIPSIS

    Traceback (most recent call last):
    ...
    RuntimeError: gamete count does not match number of diploids referring to it

Seeding a single-locus simulation from msprime
---------------------------------------------------------------------------------------------------------

Seeding a multi-locus simulation from msprime
---------------------------------------------------------------------------------------------------------

.. _fwdpp: http://molpopgen.github.io/fwdpp
