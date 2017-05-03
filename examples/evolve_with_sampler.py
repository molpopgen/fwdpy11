import fwdpy11 as fp11
from fwdpy11.model_params import SlocusParams
import fwdpy11.fitness as fp11w
import fwdpy11.wright_fisher
import numpy as np

def evolve_and_return_with_sampler(args):
    """
    Evolve function takes an initial
    population size, replicate ID number,
    RNG seed, and the type name of a recorder as
    arguments.
    """
    N,repid,seed,recorderType=args
    pop = fp11.SlocusPop(N)
    rng=fp11.GSLrng(seed)
    params=SlocusParams(
        nregions=[fp11.Region(0,1,1)],
        sregions=[fp11.ExpS(0,1,1,-0.1,1.0)],
        recregions=[fp11.Region(0,1,1)],
        gvalue=fp11w.SlocusAdditive(2.0),
        demog=np.array([N]*N,dtype=np.uint32),
        rates=(1e-3,5e-3,1e-3))
    #use the type name passed in 
    #to create an instance:
    recorder=recorderType()
    fp11.wright_fisher.evolve(rng,pop,params,recorder)
    return (repid,pop,recorder)
