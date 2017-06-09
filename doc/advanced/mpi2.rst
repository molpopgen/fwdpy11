.. _mpi2:

Parallel execution using MPI (Alternative, slower implementation)
======================================================================

Background reading:

* :ref:`mpi`

This implementation differs from that in :ref:`mpi` in that the master
process (rank 0) is coded up as a job server.  This implementation
results in a lot more interprocess communication, meaning a higher
overall system load, and slower performance.

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

    reps = [i for i in range(NREPS)]

    if rank == 0:
        conn = sqlite3.connect("output.db")
        nworkers = size - 1
        nworkers_done = 0
        ofn='output.db'
        if os.path.exists(ofn):
            os.remove(ofn)
        conn = sqlite3.connect(ofn)
        seedreps = [(i,j) for i,j in zip(seeds,reps)]

        while nworkers_done < nworkers:
            data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
            source = status.Get_source()
            tag = status.Get_tag()
            if tag == tags.READY:
                if len(seedreps) > 0:
                    comm.send(seedreps[0],dest=source,tag=tags.START)
                    seedreps.pop(0)
                    print(len(seedreps), " jobs remaining")
                else:
                    comm.send(None,dest=source,tag=tags.EXIT)
            elif tag == tags.DONE:
                print("got data from ", source)
                results = data
                df = pd.DataFrame(data)
                df.to_sql('pi',conn,if_exists='append')
            elif tag == tags.EXIT:
                nworkers_done += 1
        conn.close()
    else:
        while True:
            comm.send(None, dest=0, tag=tags.READY) 
            task = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
            tag = status.Get_tag()
            if tag == tags.START:
                N=1000
                pop = fwdpy11.SlocusPop(N)
                rng=fwdpy11.GSLrng(task[0])
                params=fwdpy11.model_params.SlocusParams(
                    nregions=[fwdpy11.Region(0,1,1)],
                    sregions=[fwdpy11.ExpS(0,1,1,-0.1,1.0)],
                    recregions=[fwdpy11.Region(0,1,1)],
                    gvalue=fwdpy11.fitness.SlocusAdditive(2.0),
                    demography=np.array([N]*10*N,dtype=np.uint32),
                    rates=(1e-3,5e-3,1e-3))

                recorder = Pi(50,task[1],rng)
                fwdpy11.wright_fisher.evolve(rng,pop,params,recorder)
                comm.send(recorder.data,dest=0,tag=tags.DONE)
            elif tag == tags.EXIT:
                break
        #Tell the master that we're done
        comm.send(None,dest=0,tag=tags.EXIT)
