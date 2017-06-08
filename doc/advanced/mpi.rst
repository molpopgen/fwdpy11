.. _mpi:

Parallel execution using MPI
======================================================================


.. code-block:: python

    from mpi4py import MPI
    import fwdpy11 as fpl1
    import fwdpy11.sampling
    import fwdpy11.model_params
    import fwdpy11.fitness
    import fwdpy11.wright_fisher
    import numpy as np

    comm = MPI.COMM_WORLD

    class Pi(object):
        """
        Calculate pi (sum of site heterozygosity)
        from NumPy matrices based on a sample.

        The result is a time series pi over time.
        """
        def __init__(self,nsam,rng):
            self.data=[]
            self.nsam=nsam
            self.rng=rng
        def __call__(self,pop):
            ind = np.random.choice(pop.N,self.nsam,replace=False)
            keys=fwdpy11.sampling.mutation_keys(pop,ind)
            if len(keys[0]) == 0: #There are no neutral variants in the sample
                self.data.append(0.)
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
                self.data.append((pop.generation,ssh))

    np.random.seed(42)

    seeds = sorted(np.random.choice(int(4e6),comm.Get_size(),replace=False))

    rank = comm.Get_rank()
    seed = seeds[rank]


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

    recorder = Pi(50,rng)
    fwdpy11.wright_fisher.evolve(rng,pop,params,recorder)

    sampler_results = comm.gather(recorder.data,root=0)

    if rank == 0:
        for i in sampler_results:
            print(len(i),i[-1])

