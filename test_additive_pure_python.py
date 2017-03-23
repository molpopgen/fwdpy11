import fwdpy11 as fp11
import fwdpy11.wright_fisher
import multiprocessing as mp
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
    def __call__(self,pop,generation):
        return
        c=Counter()
        for m in pop.mcounts:
            if m > 0:
                c[m]+=1
        #The sample of size nsam is taken by calling
        #another pybind11-based function that calls
        #fwdpp in the back-end
        sample=fp11.sample_separate(self.rng,pop,self.nsam,True);
        self.data.append((generation,c,sample))

def additive(d,g,m):
    s=0.;
    for mk in g[d.first].smutations:
        s+=m[mk].s
    for mk in g[d.second].smutations:
        s+=m[mk].s
    return max(0,1.0+s)

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
    fwdpy11.wright_fisher.evolve_regions_sampler_fitness(rng,pop,rec,additive,N,1000,0.001,0.005,0.001,nregions,sregions,recregions)
    #OMG pops are now pickle-able!!!
    return (pop,rec)

if __name__ == "__main__":
    #run 10 sims in parallel with a fitness fxn written in Python
    np.random.seed(101)
    args=[(1000,seed) for seed in np.random.randint(0,42000000,10)]
    P=mp.Pool()
    res=P.imap(evolve_and_return,args)
    P.close()
    P.join()

    for i in res:
        print(len(i[0].mutations))

