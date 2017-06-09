.. _parallel:

Running simulations in parallel
==========================================

Using Python's built-in methods for running multiple processes
-------------------------------------------------------------------------------

The general recipe is:

* Define a function to run your sims.  The function takes a `tuple` as arguments and returns a pickleable type.
* Use Python 3's concurrent.futures_ to run multiple simulations in parallel.

.. _concurrent.futures: https://docs.python.org/3/library/concurrent.futures.html

.. note::
    Be very careful not to confuse **threads** with **processes**.  The former will result in the Global Interpreter
    Lock (GIL) grinding your work to a halt.  The later spawns off work to truly inependent Python processes.  When
    reading the concurrent.futures_ documentation, focus on the sections related to processes.

.. ipython:: python
    :suppress:

    import fwdpy11
    import fwdpy11.wright_fisher
    import fwdpy11.ezparams
    import concurrent.futures
    import numpy
    import fwdpy11.fitness
    import fwdpy11.model_params
    import collections

Define an "evolve" function
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Running simulations in separate Python processes requires a function that does the work.  The arguments 
to that function (in most cases) must be a tuple of pickleable types.

Below, we define an "evolve" function that will do the following:

* Takes the population size, replicated ID numbers, random number seed, and type name of a temporal sampler/recorder as arguments.
* Initializes our population and model parameters objects
* Runs the simulation
* Returns a tuple of replicate ID, the simulated population, and an instance of the temporal sampler.

It has the following "bells and whistles":

* The temporal sampler is optional. The function relies on the fact that Python functions may take a *type name* as an argument or an *instance* of a type.  If `None` is passed in, then no recorder is used.  Else, an instance of the type name is created and used.

.. note::
    The return value from a function executed in another Python process must be pickleable.  

Here is our function:

.. ipython:: python

    def evolve(args):
        """
        Evolve function takes an initial
        population size, replicate ID number,
        RNG seed, and the type name of a recorder as
        arguments.
        """
        N,repid,seed,recorderType=args
        pop = fwdpy11.SlocusPop(N)
        rng=fwdpy11.GSLrng(seed)
        params=fwdpy11.model_params.SlocusParams(
            nregions=[fwdpy11.Region(0,1,1)],
            sregions=[fwdpy11.ExpS(0,1,1,-0.1,1.0)],
            recregions=[fwdpy11.Region(0,1,1)],
            gvalue=fwdpy11.fitness.SlocusAdditive(2.0),
            #Only simulate 10 generations so 
            #that example runs quickly:
            demography=numpy.array([N]*10,dtype=np.uint32),
            rates=(1e-3,5e-3,1e-3))
        #use the type name passed in 
        #to create an instance if applicable:
        recorder = None
        if recorderType is not None:
            recorder=recorderType()
        fwdpy11.wright_fisher.evolve(rng,pop,params,recorder)
        return (repid,pop,recorder)

Let's use this function with no recorder type.

We will use numpy to generate random number seeds for our replicates.  First, we seed numpy's rng:

.. ipython:: python

    np.random.seed(101)

Next, generate a list of arguments for our processes.  For his example, we will run ten replicates.  Our arguments for each replicate
contain four elements: population size, replicate ID, seed, and recorder type name.  As we will not use a recorder, we will pass `None`:

.. ipython:: python

    args=[(1000,repid,seed,None) for repid,seed 
        in zip(range(10),np.random.randint(0,42000000,10))]
        
In order to execute the simulations in parallel, we use a process pool with a max of 10 workers

.. ipython:: python

    with concurrent.futures.ProcessPoolExecutor(10) as pool:
        #Run our simulations and get the 
        #result back, which will be
        #the population
        for res in pool.map(evolve,args):
            print(res)

Recording the site-frequency spectrum every generation
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Let's see how to use a recorder with our evolve function.  We will define a recorder to track
the site-frequency spectrum for all mutations in the population.  The most idiomatic (and fastest
method) to do this in Python is with `collections.Counter`.:

.. ipython:: python

    class RecordSFS:
        """
        This sampler records the SFS
        for the entire pop.
        """
        def __init__(self):
            #We will record data as a list
            #of tuples:
            self.data=[]
        def __call__(self,pop):
            """
            The call operator will
            be passed in the entire
            population. You can operate 
            on it in a read-only fashion
            with no copies being made.
            Basically, you're talking
            to the C++ back-end via Python.
            """
            c=collections.Counter()
            #A population records the 
            #number of occurrences of every mutation
            #in a list called 'mcounts'
            for m in pop.mcounts:
                #mcounts can contain extinct mutations,
                #so you need to skip those. They are kept
                #around because the simulation can re-use
                #their locations in memory for new mutations.
                if m > 0:
                    c[m]+=1
            #Update our sampler's data.
            self.data.append((pop.generation,c))


Apply it as we did above.  We will store the last element of the recorder from each replicate in a dict:

.. ipython:: python

    np.random.seed(101)
    data={}
    args=[(1000,repid,seed,RecordSFS) for repid,seed 
        in zip(range(10),np.random.randint(0,42000000,10))]
    with concurrent.futures.ProcessPoolExecutor(10) as pool:
        for res in pool.map(evolve,args):
            print(res)
            data[res[0]]=res[2].data[-1]

Asynchronous execution using multiprocessing
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

You'll notice that the replicate ID numbers were output in a sorted order in the above output.  That does not have to be the case.
For example, we could use the multiprocessing_ module instead of concurrent.futures_.  The multiprocessing_ modules 'imap_unordered` function returns
results in an unpredictable order (which motivates passing a replicate ID to the evolve function).  However, the results for a given replicate will be the same
as in the previous example because we've used the same set of RNG seeds:

.. ipython:: python

    import multiprocessing as mp
    P=mp.Pool(10)
    res=P.imap_unordered(evolve,args)
    P.close()
    P.join()

Let's compare the output to what we stored from our previous set of simulations:

.. ipython:: python

    for i in res:
        print(i)
        print(i[2].data[-1] == data[i[0]])

As we can see, multiprocessing_ gives us the same results but in a random order.  Using multiprocessing_ or concurrent.futures_ is largely a matter of preference.  
The latter is Python3 only, but so is this package.  It is also more "streamlined" in its interface than multiprocessing_, which also makes it less flexible.

.. _multiprocessing: https://docs.python.org/3.5/library/multiprocessing.html
.. _concurrent.futures: https://docs.python.org/3/library/concurrent.futures.html
