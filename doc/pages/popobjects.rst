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
    # and contains our mutation
    gametes.append(fwdpy11.Gamete((2,fwdpy11.VectorUint32([]),fwdpy11.VectorUint32([0]))))

    # There is only one gamete, and so our diploid
    # will contain two copies of it:
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
