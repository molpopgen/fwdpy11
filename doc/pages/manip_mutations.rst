.. _manipulating_mutations:

Adding mutations to a population
======================================================================

To add one or more new mutations to a population, make a call to the following member functions:

* :attr:`fwdpy11.SlocusPop.add_mutations`
* :attr:`fwdpy11.MlocusPop.add_mutations`

Briefly, mutations are added to specific gametes of individuals.  For example, to add a new mutation
to the second gamete of the third individual in the population:

.. ipython:: python

    import fwdpy11

    pop = fwdpy11.SlocusPop(1000)
    new_mutations = fwdpy11.VecMutation()
    new_mutations.append(fwdpy11.Mutation(1.,0,0,0,0))
    print(new_mutations[0])
    new_mutation_keys = pop.add_mutations(new_mutations, [3], [1])
    assert(len(new_mutation_keys)==1)
    assert(pop.mutations[new_mutation_keys[0]].neutral is True)
    assert(pop.mcounts[new_mutation_keys[0]] == 1)

See the unit test file tests/test_add_mutations.py for more examples.

.. note::

    Attempting to add a mutation at an already-mutation position will raise an exception,
    because you are violating the infinitely-many sites assumption of the simulation engine.
    To guard against this, you may check the current set of mutation positions via
    :attr:`fwdpy11.Population.mut_lookup`.  See :ref:`mpos` for details.

Changing the effect size of a mutation
======================================================================

To change the effect size of a mutation, call :func:`fwdpy11.util.change_effect_size`:

.. ipython:: python

    import fwdpy11.util

    # Make our mutation non-neutral, via
    # changing it to have a vector of nonzero effect sizes:
    fwdpy11.util.change_effect_size(pop, new_mutation_keys[0], new_esizes = [-1.,2.])
    print(pop.mutations[new_mutation_keys[0]].neutral)
    print(pop.mutations[new_mutation_keys[0]].esizes)

    # Change it back to a mutation w/constant effect size:
    fwdpy11.util.change_effect_size(pop, new_mutation_keys[0])
    print(pop.mutations[new_mutation_keys[0]].neutral)
    print(pop.mutations[new_mutation_keys[0]].esizes)

The above example makes heavy use of default parameters values in :func:`fwdpy11.util.change_effect_size`.  See that
function's docstring for details.
