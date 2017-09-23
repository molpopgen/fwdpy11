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

.. testcode::

    import fwdpy11

    mutations = [fwdpy11.Mutation(0.1,-0.01,1.0,0,0)]
    gametes = [fwdpy11.Gamete((2,fwdpy11.VectorUint32([]),fwdpy11.VectorUint32([0])))]
    diploids = [fwdpy11.SingleLocusDiploid(0,0)]
    pop = fwdpy11.SlocusPop(diploids, gametes, mutations)
