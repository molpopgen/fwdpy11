import fwdpy11 as fp11
import fwdpy11.fitness as fp11w
import fwdpy11.wright_fisher
import numpy as np

def evolve_and_return_with_sampler(args):
    print("here")
    N,repid,seed,recorderType=args
    pop = fp11.Spop(N)
    rng=fp11.GSLrng(seed)
    nregions=[fp11.Region(0,1,1)]
    sregions=[fp11.ExpS(0,1,1,-0.1,1.0)]
    recregions=nregions
    fitness=fp11w.SpopAdditive(2.0)
    nlist=np.array([N]*N,dtype=np.uint32)
    recorder=recorderType()
    fp11.wright_fisher.evolve_regions_sampler_fitness(rng,pop,
            nlist,0.001,0.005,0.001,nregions,sregions,recregions,fitness,recorder)
    return (repid,pop,recorder)
