.. _pickling_pops:

Pickling populations
==========================================

Population objects may be pickled, enabling running multiple replicates using multiprocessing_ or concurrent.futures_.
Other applications include saving the state of a population in memory in order to restart evolution from the same
starting conditions.

The pickling is implemented via fwdpp's binary serialization methods.  This is faster and more space-efficient than
natively pickling the Python representation of the underlying C++ types.  

In order to pickle a file, you must use the latest pickling protocol.

Some general notes:

* We are trying hard to keep the binary formats constant, and to use "magic" numbers to read in older format versions
  whenever possible.  However, in 0.1.5, we had to make a clean break, and attempting to unpickle populations from
  previous versions will result in an exception.

.. note::
    If you are pickling populations to a file, you should be aware that you can (and probably should) compress the output.
    We have found that lzma_ compression results in the smallest files, but takes the longest to compress.  The gzip_
    module is a good balance between speed and final file size.

.. warning::
    If you are simulating very large populations and/or a very large number of mutations, it is possible to 
    run out of memory during pickling.  When this occurs, exceptions will be thrown from fwdpp and translated
    to a Python RuntimeError.  The cases where we have experienced such failures is when
    simulating on the order of 10 megabases of variation under the Tennessen_ model of European demography. The 
    last few time points of that model involve very large population sizes.
    
    If you wish to serialize a population to file, an alternative to pickling is described in :ref:`binary_pops`.

In fwdpy11 0.3.0, :class:`fwdpy11.SlocusPop` has two new functions that avoid high memory consumption while pickling.
The cost is that pickling is slower.  The functions are:

* :func:`fwdpy11.SlocusPop.pickle_to_file`, which is an instance method
* :func:`fwdpy11.SlocusPop.load_from_pickle_file`, which is a static method

Please notes that these functions must be used in a coordinated manner!  See documentation for details.

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
    params = fp11.model_params.ModelParams(**p)
    rng = fp11.GSLrng(42)
    wf.evolve(rng,pop,params)
    #Pickle the pop in memory.
    #The -1 means use latest protocol
    ppop = pickle.dumps(pop,-1)
    #Unpickle to create a new pop:
    pop2 = pickle.loads(ppop)
    print(pop==pop2)

.. testoutput::

    True

.. _multiprocessing: https://docs.python.org/3/library/multiprocessing.html
.. _concurrent.futures: https://docs.python.org/3/library/concurrent.futures.html
.. _lzma: https://docs.python.org/3/library/lzma.html
.. _gzip: https://docs.python.org/3/library/gzip.html
.. _Tennessen: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3708544/
