.. _pickling_pops:

Pickling populations
==========================================

Population objects may be pickled, enabling running multiple replicates using multiprocessing_ or concurrent.futures_.
Other applications include saving the state of a population in memory in order to restart evolution from the same
starting conditions.

The pickling is implemented via fwdpp's binary serialization methods.  This is faster and more space-efficient than
natively pickling the Python representation of the underlying C++ types.  

Full support for pickling is only possible with Python 3.  Under Python 2, you may pickle/unpickle to/from memory just
fine, but you cannot unpickle from a file.

In order to pickle a file, you must also use the latest pickling protocol.

.. note::
    If you are pickling populations to a file, you should be aware that you can (and probably should) compress the output.
    We have found that lzma_ compression is ideal when pickling fwdpy11 population objects.

.. warning::
    If you are simulating very large populations and/or a very large number of mutations, it is possible to 
    run out of memory during pickling.  When this occurs, exceptions will be thrown from fwdpp and translated
    to a Python RuntimeError.  The cases where we have experienced such failures is when
    simulating on the order of 10 megabases of variation under the Tennessen_ model of European demography. The 
    last few time points of that model involve very large population sizes.
    
    Future versions of fwdpy11 may allow direct serialization to a file, which 
    avoids this problem at the cost of a less convenient API. A workaround is to pickle a sample from the 
    population, rather than the whole thing.  In this case, you may also wish to pickle the fixations, etc.,
    or whatever additional data you may need.  

.. testcode::

    import fwdpy11 as fp11
    import fwdpy11.wright_fisher as wf
    import fwdpy11.ezparams
    import fwdpy11.model_params
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
.. _lzma: https://docs.python.org/3/library/lzma.html
.. _Tennessen: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3708544/
