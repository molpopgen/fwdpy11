.. _parentage:

Tracking parentage
======================================================================

.. versionadded:: 0.1.4

This example shows one way to track parentage during a simulation.  We'll run a single-deme simulation and use
:attr:`fwdpy11.fwdpy11_types.SlocusPop.popdata_user` to record the complete pedigree for the population over time.

The mechanics are quite simple.  We define a custom recorder that gathers the parent data
(:attr:`fwdpy11.fwdpp_types.SingleLocusDiploid.parental_data`) into a list each generation.  

A real-world example would have to do more work than this.  You would probably need to assign unique individual IDs each
generation, etc., and get the data into some standard format (PLINK's ped/fam, for example), in order to do interesting
work with the pedigree.

.. ipython:: python

    import fwdpy11
    import fwdpy11.fitness
    import fwdpy11.model_params
    import fwdpy11.ezparams
    import fwdpy11.wright_fisher
    import numpy as np

    class RecordParents(object):
        def __call__(self,pop):
            parents = [i.parental_data for i in pop.diploids]
            pop.popdata_user.append(parents)

    rng = fwdpy11.GSLrng(42)

    theta,rho = 100.0,100.0
    pop = fwdpy11.SlocusPop(10)
    pop.popdata_user = []
    pdict = fwdpy11.ezparams.mslike(pop,simlen=pop.N)
    params = fwdpy11.model_params.SlocusParams(**pdict)

    r = RecordParents()
    fwdpy11.wright_fisher.evolve(rng,pop,params,r)
    for i in pop.popdata_user:
        print(i)
