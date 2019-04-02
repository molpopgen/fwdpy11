Overview of types
=============================================

This section is a conceptual overview of how a population of diploids is represented
by :class:`fwdpy11.DiploidPopulation`.

.. note::

    :class:`fwdpy11.DiploidPopulation` inherits from the base class
    :class:`fwdpy11.Population`, and thus inherits its public attributes.
    Below, you will see attribute references to the base class. Keep
    in mind that those attributes are inherited.

* A diploid is made up of two haploid genomes, one from each parent.
* A haploid genome is represented by :class:`fwdpy11.HaploidGenome`.
* Haploid genomes are stored in :attr:`fwdpy11.Population.haploid_genomes`.
* A diploid is represented by :class:`fwdpy11.DiploidGenotype`, whose attributes
  :attr:`fwdpy11.DiploidGenotype.first` and :attr:`fwdpy11.DiploidGenotype.second`
  refer to *indexes* in :attr:`fwdpy11.Population.haploid_genomes`
* Instances of :class:`fwdpy11.DiploidGenotype` are stored in :attr:`fwdpy11.DiploidPopulation.diploids`

To go from haploid genomes to their mutations:

* Mutations are represented by :class:`fwdpy11.Mutation`
* Mutation instances are stored in :attr:`fwdpy11.Population.mutations`.
* Instances of haploid genomes store *indexes* into :attr:`fwdpy11.Population.mutations`.
* These indexes are stored *separately* for neutral and selected mutations in
  :attr:`fwdpy11.HaploidGenome.mutations` and :attr:`fwdpy11.HaploidGenome.smutations`,
  respectively.

.. note::

    The attributes :attr:`fwdpy11.Population.haploid_genomes`, :attr:`fwdpy11.Population.mutations`,
    and :attr:`fwdpy11.DiploidPopulation.diploids` mostly behave as regular Python lists.  However,
    they are actually C++ containers and some magic has been done to allow you to access their
    data very efficiently.

Let's take a look at the population simulated in :ref:`introexample`:

.. ipython:: python

    print(pop.diploids[0].first, pop.diploids[0].second)

.. ipython:: python

    for i in (pop.diploids[0].first, pop.diploids[0].second):
        print(pop.haploid_genomes[i].smutations,
              type(pop.haploid_genomes[i].smutations))
              pop.haploid_genomes[i].smutations.dtype)

So we see that mutation indexes are stored in numpy arrays.

.. note::

    :attr:`fwdpy11.HaploidGenome.mutations` is empty in simulations 
    with tree sequences!  Neutral variants are added after-the-fact
    and are processed entirely from the tree sequence.

Let's take a look at the mutations for one of the genomes:

.. ipython:: python

    for k in pop.haploid_genomes[pop.diploids[0].first].smutations:
        print(pop.mutations[k].pos, pop.mutations[k].g, pop.mutations[k].s)

