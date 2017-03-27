import fwdpy11 as fp11
import fwdpy11.wright_fisher
import multiprocessing as mp
from collections import Counter
import numpy as np

class RecordSFS:
    def __init__(self):
        self.data=[]
    def __call__(self,pop):
        c=Counter()
        for m in pop.mcounts:
            if m > 0:
                c[m]+=1
        self.data.append((pop.generation,c))

def evolve_and_return(args):
    N,seed=args
    pop = fp11.Spop(N)
    rec=RecordSFS()
    rng=fp11.GSLrng(seed)
    fwdpy11.wright_fisher.evolve(pop,rng,1000,10000,0.001,0.001)
    #OMG pops are now pickle-able!!!
    print(len(pop.mutations),len(pop.mcounts))
    return (pop,rec)

if __name__ == "__main__":
    args=[(1000,seed) for seed in np.random.randint(0,42000000,10)]
    print(args)
    P=mp.Pool()
    res=P.imap(evolve_and_return,args)
    P.close()
    P.join()
