import fwdpy11 as fp11
import fwdpy11.fitness as fp11w
import fwdpy11.sampling as fp11s
import fwdpy11.wright_fisher
import concurrent.futures
from collections import Counter
import numpy as np
import math

class RecordSFSandSample:
    """
    This sampler records the SFS
    for the entire pop and a sample
    of size nsam.
    """
    def __init__(self,rng,nsam):
        self.data=[]
        self.rng=rng
        self.nsam=nsam
    def __call__(self,pop):
        c=Counter()
        for m in pop.mcounts:
            if m > 0:
                c[m]+=1
        #The sample of size nsam is taken by calling
        #another pybind11-based function that calls
        #fwdpp in the back-end
        sample=fp11s.sample_separate(self.rng,pop,self.nsam,True);
        self.data.append((pop.generation,c,sample))

def evolve_and_return(args):
    """
    Caveat: C++-based fitness types
    are not pickleable, and thus cannot
    be sent in this way. Have to think of a workaround
    """
    N,seed=args
    pop = fp11.Spop(N)
    rng=fp11.GSLrng(seed)
    rec=RecordSFSandSample(rng,10)
    nregions=[fp11.Region(0,1,1)]
    sregions=[fp11.ExpS(0,1,1,-0.1,1.0)]
    recregions=nregions
    fitness=fp11w.SpopAdditive(2.0)
    nlist=np.array([N]*N,dtype=np.uint32)
    fwdpy11.wright_fisher.evolve_regions_sampler_fitness(rng,pop,
            nlist,0.001,0.005,0.001,nregions,sregions,recregions,fitness,rec)
    return (pop,rec)

if __name__ == "__main__":
    #run 10 sims in parallel with a fitness fxn written in Python
    np.random.seed(101)
    args=[(1000,seed) for seed in np.random.randint(0,42000000,10)]

    with concurrent.futures.ProcessPoolExecutor(10) as pool:
        for res in pool.map(evolve_and_return,args):
            print (type(res))
            print(res[0].mcounts,len(res[1].data[0]))
            print(res[1].data[-1])
            print (res[0].gametes)
            for dip in res[0].diploids:
                for mk in res[0].gametes[dip.first].smutations:
                    print( res[0].mutations[mk])
