Running simulations in parallel
==========================================

.. code-block:: python

    import fwdpy11 as fp11
    import fwdpy11.wright_fisher as wf
    #concurrent.futures is Python 3 only
    import concurrent.futures as cf
    import numpy as np

    def evolve_and_return(args):
        """
        This function runs our simulation.
        The input arguments come in a tuple,
        which is required by many of Python's
        functions for execution in separate processes.

        For this function, the arguments are the population
        size and a random number seed.
        """
        N,seed=args
        #Construct as single-deme object
        #with N diploids
        pop = fp11.Spop(N)
        #Initialize a random number generator
        rng=fp11.GSLrng(seed)
        sregions=[fp11.ExpS(0,1,1,-0.1,1.0)]
        #Evolve the pop for 10N generations
        #and several default parameters will be
        #used.
        fp11.wright_fisher.evolve(rng,pop,sregions=sregions)
        #The population is picklable, and so
        #we can return it from another process
        return pop 

    if __name__ == "__main__":
        #init global rng seed
        np.random.seed(101)
        #Generate a list of arguments for our processes.
        #We generation 10 random number seeds
        args=[(1000,seed) for seed in np.random.randint(0,42000000,10)]
        #Use a process pool with a max of 10 workers
        with cf.ProcessPoolExecutor(10) as pool:
            #Run our simulations and get the 
            #result back, which will be
            #the population
            for res in pool.map(evolve_and_return,args):
                print(type(res))
                
