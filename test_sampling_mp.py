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
        c=Counter()
        for m in pop.mcounts:
            if m > 0:
                c[m]+=1
        #The sample of size nsam is taken by calling
        #another pybind11-based function that calls
        #fwdpp in the back-end
        sample=fp11.sample_separate(self.rng,pop,self.nsam,True);
        self.data.append((generation,c,sample))

#two custom fitness functions:

def neutral_fitness(a,b,c):
    return 1

def exp_decline_with_ttl(dip,gametes,mutations):
    ttl=0.0
    for mk in gametes[dip.first].smutations:
        ttl += mutations[mk].s
    for mk in gametes[dip.second].smutations:
        ttl += mutations[mk].s
    return math.exp(-(math.fabs(ttl)))

def evolve_and_return(args):
    """
    Caveat: C++-based fitness types
    are not pickleable, and thus cannot
    be sent in this way. Have to think of a workaround
    """
    N,seed,fitness=args
    pop = fp11.Spop(N)
    rng=fp11.GSLrng(seed)
    rec=RecordSFSandSample(rng,10)
    nregions=[fp11.Region(0,1,1)]
    sregions=[fp11.ExpS(0,1,1,-0.1,1.0)]
    recregions=nregions
    fitness=fp11.SpopAdditive(2.0).callback
    fwdpy11.wright_fisher.evolve_regions_sampler_fitness(rng,pop,
            N,1000,0.001,0.005,0.001,nregions,sregions,recregions,fitness,rec)
    print(len(pop.mutations))
    #OMG pops are now pickle-able!!!
    #return (pop,rec)

if __name__ == "__main__":
    #run 10 sims in parallel with a fitness fxn written in Python
    np.random.seed(101)
    args=[(1000,seed,exp_decline_with_ttl) for seed in np.random.randint(0,42000000,10)]
    evolve_and_return(args[0])
    #P=mp.Pool()
    #res=P.imap(evolve_and_return,args)
    #P.close()
    #P.join()

    #for r in res:
    #    print(r[0].mcounts,len(r[1].data[0]))
    #    print(r[1].data[-1])
    #    print (r[0].gametes)
    #    for dip in r[0].diploids:
    #        for mk in r[0].gametes[dip.first].mutations:
    #            print( r[0].mutations[mk])
