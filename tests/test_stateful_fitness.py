#We use cppimport to build the module.
#For the sake of unit testing, we force
#a rebuild every time, but that clearly
#wouldn't be needed for normal use.
import cppimport
cppimport.force_rebuild()
cppimport.set_quiet(False)
snowdrift = cppimport.imp("snowdrift")
import unittest
import numpy as np
import fwdpy11 as fp11
import fwdpy11.temporal_samplers as fp11ts
import fwdpy11.wright_fisher
def evolve_snowdrift(args):
    """
    We write the function taking a tuple 
    out of habit, simplifying later
    integration with multiprocessing or
    concurrent.futures.
    """
    N,seed=args
    #Construct as single-deme object
    #with N diploids
    pop = fp11.Spop(N)
    #Initialize a random number generator
    rng=fp11.GSLrng(seed)
    sregions=[fp11.ExpS(0,1,1,-0.1,1.0)]
    recregions=[fp11.Region(0,1,1)]
    fitness = snowdrift.SpopSnowdrift(0.2,-0.2,1,-2)
    #evolve for 100 generations so that unit tests are
    #fast
    nlist = np.array([N]*100,dtype=np.uint32)
    recorder = fp11ts.RecordNothing()
    fp11.wright_fisher.evolve_regions_sampler_fitness(rng,pop,nlist,0.0,0.0025,0.001,
            [],sregions,recregions,fitness,recorder,prune_all_fixations=False)
    #return our pop
    return pop

class testSnowdrift(unittest.TestCase):
    def test_create(self):
        f = snowdrift.SpopSnowdrift(1,-1,0.1,0.2);
    def test_evolve(self):
        p = evolve_snowdrift((1000,42))

if __name__=="__main__":
    unittest.main()
