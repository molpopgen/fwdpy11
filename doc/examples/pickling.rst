Pickling populations
==========================================

Population objects may be pickled, enabling running multiple replicates using multiprocessing_ or concurrent.futures_.
Other applications include saving the state of a population in memory in order to restart evolution from the same
starting conditions.

The pickling is implemented via fwdpp's binary serialization methods.  This is faster and more space-efficient than
natively pickling the Python representation of the underlying C++ types.  

Full support for pickling is only possible with Python 3.  Under Python 2, you may pickle/unpickle to/from memory just
fine, but you cannot unpickle from a file.

In order to picked a file, you must also use the latest pickling protocol.

.. testcode::

    import fwdpy11 as fp11
    import fwdpy11.wright_fisher as wf
    import fwdpy11.ezparams
    import numpy as np
    import pickle

    pop = fp11.SlocusPop(1000)
    p = fwdpy11.ezparams.mslike(pop,
        simlen=100,
        rates=(1e-3,0.,1e-3))
    params = fp11.model_params.SlocusParams(**p)
    rng = fp11.GSLrng(42)
    wf.evolve(rng,pop,params)
    #Pickle the pop in memory.
    #The -1 means use latest protocol
    ppop = pickle.dumps(pop,-1)
    #Unpickle to create a new pop:
    pop2 = pickle.loads(ppop)
    print(type(pop2))
    print(pop==pop2)

.. testoutput::

    <class 'fwdpy11.fwdpy11_types.SlocusPop'>
    True

.. _multiprocessing: https://docs.python.org/3/library/multiprocessing.html
.. _concurrent.futures: https://docs.python.org/3/library/concurrent.futures.html
