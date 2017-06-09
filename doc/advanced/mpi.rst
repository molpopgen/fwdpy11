.. _mpi:

Parallel execution using MPI
======================================================================

Background reading:

* :ref:`parallel`

An alternative way to run replicate simulations in parallel is to use the Message Passing Interface, or MPI.
MPI allows different processes to communicate back and forth with one another.  This communication can happen
within and between different computers, meaning that you can farm replicates out to multiple compute nodes on a
cluster.

In Python, MPI is available via the mpi4py_ module.  

.. note:: 
    In general, MPI is a complex beast, and it is beyone the scope
    of this manual do teach it.  Further, its installation/configuration can be a bit, well, sensitive, for lack of a better
    term. On your own machine, it will most likely "just work" when installed via Anaconda.  However, 
    on a managed compute cluster, you may need the assistance of your support folks.

Similar to our previous example (:ref:`parallel`) based on concurrent.futures_, a big advantage of running replicates
like this is the ability to collect all results into a common output file.

One common MPI idiom is to have one process be devoted to farming out work and receiving output from other processes.
Here, we will use this "master" task to collect time series data from simulations.

One critical issue controlling MPI performance is the extent of interprocess communication.  In the example below, we minimize the communication by having each task with rank > 0 figure out which replicates it will run and which seeds it will use.  We accomplish this calculation using Python "fancy indexing" within each task with rank > 0. 

.. note::
    An alternative to having each task calculate its set of operations is to have the task with rank 0
    act like a "job server", meaning that it asks all other tasks if they are ready to work, sends them data
    if they are, gets the data back when they're done, and then tells them to exit.  Ony my 40 core 
    machine, such an implementation takes twice as long to execute as the example below.  For the sake of 
    completeness, though, you can see it in :ref:`mpi2`

Our implementation is such that the task with rank 0 is the "master", and all other tasks (rank > 0) run the simulations.  Thus, your MPI world must have a size of at least 2:

.. code-block:: bash

    #Pass -n >= 2:
    mpirun -n 2 python3 my_fwdpy11_mpi_script.py

The example below does the following:

1. 1,000 simulations are run.
2. Each simulation records :math:`\theta_\pi` in a sample of 50 diploids each generation.
3. The time series for each replicate is returned to the master process.
4. The master process converts all data into a Pandas DataFrame and dumps the results to an sqlite3 database. 

.. code-block:: python

    #implementation based on
    #https://github.com/jbornschein/mpi4py-examples/blob/master/09-task-pull.py

    from mpi4py import MPI
    import fwdpy11 as fpl1
    import fwdpy11.sampling
    import fwdpy11.model_params
    import fwdpy11.fitness
    import fwdpy11.wright_fisher
    import numpy as np
    import pandas as pd
    import sqlite3,os
    import time

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    status = MPI.Status()

    #Set the global RNG seed
    np.random.seed(42)

    #We'll simulate 20 replicate
    #populations
    NREPS=1000

    def enum(*sequential, **named):
        """Handy way to fake an enumerated type in Python
        http://stackoverflow.com/questions/36932/how-can-i-represent-an-enum-in-python
        """
        enums = dict(zip(sequential, range(len(sequential))), **named)
        return type('Enum', (), enums)

    # Define MPI message tags
    tags = enum('READY', 'DONE', 'EXIT', 'START')

    def empty_array():
        return np.array([],
                dtype=[('repid',np.int32),
                    ('generation',np.int32),
                    ('pi',np.float)])

    class Pi(object):
        """
        Calculate pi (sum of site heterozygosity)
        from NumPy matrices based on a sample.

        The result is a time series pi over time.
        """
        def __init__(self,nsam,repid,rng):
            self.data=empty_array()
            self.repid=repid
            self.nsam=nsam
            self.rng=rng
        def __call__(self,pop):
            ind = np.random.choice(pop.N,self.nsam,replace=False)
            keys=fwdpy11.sampling.mutation_keys(pop,ind)
            if len(keys[0]) == 0: #There are no neutral variants in the sample
                self.data=np.concatenate([self.data,
                    np.array([(self.repid,pop.generation,0.0)],
                        dtype=self.data.dtype)])
            else:
                neutral_sorted_keys=[i for i in sorted(keys[0],
                    key=lambda x,m=pop.mutations: m[x[0]].pos)]
                dm = fwdpy11.sampling.haplotype_matrix(pop,ind,
                    neutral_sorted_keys,keys[1])
                n = dm.neutral()
                n=n.reshape(dm.nrow,len(neutral_sorted_keys))
                colsums=n.sum(0)
                ssh = 0.
                for i in colsums:
                    ssh += i*(float(2*self.nsam)-i)
                ssh *= 2.
                ssh /= float((2*self.nsam)*(2*self.nsam-1))
                self.data=np.concatenate([self.data,
                    np.array([(self.repid,pop.generation,ssh)],
                        dtype=self.data.dtype)])


    #Generate seeds for each replicate
    seeds = np.random.choice(int(4e6),
        NREPS,
        replace=False)

    #In general, we want to simulate more
    #replicates than the number of cores 
    #available to the MPI system. We'll
    #use fancy indexing to split up our
    #seeds so that each process gets 
    #an approximately even amount of work 
    #to do.
    seeds_for_ranks = [seeds[i::comm.Get_size()-1] for i in range(comm.Get_size()-1)] 
    #Do the same trick to assign a 
    #"replicate ID number" to each task:
    reps = [i for i in range(NREPS)]
    reps_for_ranks =[reps[i::comm.Get_size()-1] for i in range(comm.Get_size()-1)] 

    if rank == 0:
        conn = sqlite3.connect("output.db")
        nworkers = size - 1
        nworkers_done = 0
        ofn='output.db'
        if os.path.exists(ofn):
            os.remove(ofn)
        conn = sqlite3.connect(ofn)
        while nworkers_done < nworkers:
            data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
            source = status.Get_source()
            tag = status.Get_tag()
            if tag == tags.DONE:
                results = data
                df = pd.DataFrame(data)
                df.to_sql('pi',conn,if_exists='append')
            elif tag == tags.EXIT:
                nworkers_done += 1
        conn.close()
    else:
        #These are the seeds and 
        #repids for this specific task
        seeds = seeds_for_ranks[rank-1]
        repids = reps_for_ranks[rank-1] 
        #For each replicate, run a simulation:
        for seed,repid in zip(seeds,repids):
            N=1000
            pop = fwdpy11.SlocusPop(N)
            rng=fwdpy11.GSLrng(seed)
            params=fwdpy11.model_params.SlocusParams(
                nregions=[fwdpy11.Region(0,1,1)],
                sregions=[fwdpy11.ExpS(0,1,1,-0.1,1.0)],
                recregions=[fwdpy11.Region(0,1,1)],
                gvalue=fwdpy11.fitness.SlocusAdditive(2.0),
                demography=np.array([N]*10*N,dtype=np.uint32),
                rates=(1e-3,5e-3,1e-3))

            recorder = Pi(50,repid,rng)
            fwdpy11.wright_fisher.evolve(rng,pop,params,recorder)
            comm.send(recorder.data,dest=0,tag=tags.DONE)

        #Tell the master that we're done
        comm.send(None,dest=0,tag=tags.EXIT)

.. _concurrent.futures: https://docs.python.org/3/library/concurrent.futures.html
.. _mpi4py: http://mpi4py.readthedocs.io
