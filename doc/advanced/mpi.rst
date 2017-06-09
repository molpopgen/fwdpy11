.. _mpi:

Parallel execution using MPI
======================================================================

The example below does the following:

1. 1,000 simulations are run
2. Each simulation records :math:`\theta_\pi` in a sample of 50 dipliods each generation.
3. The time series for each replicate is pickled to a file 
4. Once all replicates are done, the MPI task with rank 0 collates all the pickled data into a single sqlite3 database.

.. code-block:: python

    from mpi4py import MPI
    import fwdpy11 as fpl1
    import fwdpy11.sampling
    import fwdpy11.model_params
    import fwdpy11.fitness
    import fwdpy11.wright_fisher
    import numpy as np
    import pandas as pd
    import sqlite3
    import gzip,pickle,os

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    #We'll simulate 1000 replicate
    #populations
    NREPS=1000

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

    #Set the global RNG seed
    np.random.seed(42)


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
    seeds_for_ranks = [seeds[i::comm.Get_size()] for i in range(comm.Get_size())] 
    #Do the same trick to assign a 
    #"replicate ID number" to each task:
    reps = [i for i in range(NREPS)]
    reps_for_ranks =[reps[i::comm.Get_size()] for i in range(comm.Get_size())] 

    #These are the seeds and repids 
    #for this specific rank
    seeds = seeds_for_ranks[rank]
    repids = reps_for_ranks[rank] 

    #Create output file for this rank
    ofn="rank" + str(rank) + ".gz"
    if os.path.exists(ofn):
        os.remove(ofn)
    of = gzip.open(ofn,'ab')
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
        pickle.dump(recorder.data,of)   
        pop.clear()
        pop=None
        del pop
        recorder.data = None
        recorder = None
        del recorder
    of.close()

    #Return all file names to rank-0 process:
    FILES = comm.gather(ofn,root=0)

    if rank == 0:
        dbname='output.db'
        if os.path.exists(dbname):
            os.remove(dbname)
        conn = sqlite3.connect(dbname)
        #Go through all .gz files, unpickle
        #data, and collect it all in
        #and sqlite3 database:
        for fi in FILES:
            with gzip.open(fi,"rb") as f:
                while True:
                    try:
                        x=pickle.load(f)
                        df=pd.DataFrame(x)
                        df.to_sql('pi',conn,if_exists='append')
                    except:
                        break
                #clean up our temp files:
                os.remove(fi)
        conn.close()
