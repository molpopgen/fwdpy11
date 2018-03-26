.. _parentage:

Tracking parentage
======================================================================

.. versionadded:: 0.1.4

Background reading:

* :ref:`recorders`

This example shows one way to track parentage during a simulation.  We'll run a single-deme simulation and use
:attr:`fwdpy11.SlocusPop.popdata_user` to record the complete pedigree for the population over time.

The mechanics are quite simple.  We define a custom recorder that gathers the parent data
(:attr:`fwdpy11.SingleLocusDiploid.parental_data`) into a list each generation.  We will couple the parent
data along with the individual's label field.  This coupling facilitates later post-processing.

For convenience, we coerce the data into a `collections.namedtuple`, which allows final storage into a
`pandas.DataFrame`.

A real-world example would have to do more work than this.  You would probably need to assign unique individual IDs each
generation, etc., and get the data into some standard format (PLINK's ped/fam, for example), in order to do interesting
work with the pedigree.  Additionally, you probably want to thin the data down to remove individuals that left no
offspring, etc. 

.. ipython:: python

    import fwdpy11
    import fwdpy11.fitness
    import fwdpy11.model_params
    import fwdpy11.ezparams
    import fwdpy11.wright_fisher
    from collections import namedtuple
    import pandas as pd

    PedigreeData = namedtuple('PedigreeData',['gen','id','p1','p2'])

    # Our recorder is very simple:
    class RecordParents(object):
        def __call__(self,pop):
            parents = [PedigreeData(pop.generation,i.label,*i.parental_data) for i in pop.diploids]
            pop.popdata_user.extend(parents)

    rng = fwdpy11.GSLrng(42)

    pop = fwdpy11.SlocusPop(10)

    # It is important to initialize
    # popdata_user, whose default
    # value is None:
    pop.popdata_user = []
    pdict = fwdpy11.ezparams.mslike(pop,simlen=pop.N)
    params = fwdpy11.model_params.SlocusParams(**pdict)

    r = RecordParents()
    fwdpy11.wright_fisher.evolve(rng,pop,params,r)

    # Store the data as a DataFrame
    df = pd.DataFrame.from_records(pop.popdata_user,columns=PedigreeData._fields)
    print(df.head())
