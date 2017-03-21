import fwdpy11 as fp11
import fwdpy11.wright_fisher
import multiprocessing as mp
from collections import Counter
import numpy as np

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

def evolve_and_return(args):
    N,seed=args
    pop = fp11.Spop(N)
    rng=fp11.GSLrng(seed)
    rec=RecordSFSandSample(rng,10)
    fwdpy11.wright_fisher.evolve(pop,rng,1000,10000,0.001,0.001,rec)
    #OMG pops are now pickle-able!!!
    return (pop,rec)

if __name__ == "__main__":
    args=[(1000,seed) for seed in np.random.randint(0,42000000,10)]
    print(args)
    P=mp.Pool()
    res=P.imap(evolve_and_return,args)
    P.close()
    P.join()

    for r in res:
        print(r[0].mcounts,len(r[1].data[0]))
        print(r[1].data[-1])
        print (r[0].gametes)
        for dip in r[0].diploids:
            for mk in r[0].gametes[dip.first].mutations:
                print( r[0].mutations[mk])
