import fwdpy11 as fp11
import fwdpy11.wright_fisher
import multiprocessing as mp
from collections import Counter
import numpy as np
import math,pickle

class PopulationPickleJar:
    """
    This sampler records the SFS
    for the entire pop and a sample
    of size nsam.
    """
    def __init__(self,fn):
        self.fn=fn
    def __call__(self,pop,generation):
        with open(self.fn,"ab") as f:
            pickle.dump((generation,pop),f,-1)
            

def evolve_and_return(args):
    """
    Caveat: C++-based fitness types
    are not pickleable, and thus cannot
    be sent in this way. Have to think of a workaround
    """
    N,seed,ofn=args
    pop = fp11.Spop(N)
    rng=fp11.GSLrng(seed)
    rec=PopulationPickleJar(ofn)
    nregions=[fp11.Region(0,1,1)]
    sregions=[fp11.ExpS(0,1,1,-0.1,1.0)]
    recregions=nregions
    fitness=fp11.SpopAdditive(2.0)
    #Fool pybind11 into seeing
    #existing references
    fwdpy11.wright_fisher.evolve_regions_sampler_fitness(rng,pop,
            N,1000,0.001,0.005,0.001,nregions,sregions,recregions,fitness,rec)
    print(len(pop.mutations))
    #OMG pops are now pickle-able!!!
    #return (pop,rec)

if __name__ == "__main__":
    #run 10 sims in parallel with a fitness fxn written in Python
    np.random.seed(101)
    args=[(1000,seed,"test.pickle") for seed in np.random.randint(0,42000000,10)]
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
