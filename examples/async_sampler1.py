import fwdpy11 as fp11
import fwdpy11.fitness as fp11w
import fwdpy11.wright_fisher
import multiprocessing as mp
from collections import Counter
import numpy as np
from evolve_with_sampler import evolve_and_return_with_sampler
import unittest

class RecordSFS:
    """
    This sampler records the SFS
    for the entire pop.
    """
    def __init__(self):
        self.data=[]
    def __call__(self,pop):
        """
        The call operator will
        be passed in the entire
        population. You can operate 
        on it in a read-only fashion
        with no copied being made.
        Basically, you're talking
        to the C++ back-end via Python.
        """
        #Use Python's collections.Counter
        #to record the SFS
        c=Counter()
        #A population records the 
        #number of occurrences of every mutation
        #in a list called 'mcounts'
        for m in pop.mcounts:
            #mcounts can contain extinct mutations,
            #so you need to skip those. They are kept
            #around because the simulation can re-use
            #their locations in memory for new mutations.
            if m > 0:
                c[m]+=1
        #Update our sampler's data.
        self.data.append((pop.generation,c))


if __name__ == "__main__":
    np.random.seed(101)
    args=[(1000,repid,seed,RecordSFS) for
            repid,seed in
            zip(range(10),np.random.randint(0,42000000,10))]
    print(args)
    P=mp.Pool()
    res=P.imap_unordered(evolve_and_return_with_sampler,args)
    P.close()
    P.join()

    for i in res:
        print(i)
