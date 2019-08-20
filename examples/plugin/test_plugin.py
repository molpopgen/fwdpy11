import fwdpy11
import numpy as np
import gvalue_recorder

N = 1000

pdict = {'demography': np.array([N]*100),
         'nregions': [],
         'sregions': [fwdpy11.GaussianS(0, 1, 1, 0.25)],
         'recregions': [fwdpy11.Region(0, 1, 1)],
         'rates': (0., 5e-3, 5e-3),
         'gvalue': fwdpy11.Additive(2., fwdpy11.GSS(VS=1, opt=1)),
         'prune_selected': False
         }

params = fwdpy11.ModelParams(**pdict)

pop = fwdpy11.DiploidPopulation(N, 1.0)

rng = fwdpy11.GSLrng(42)

gvalues = []
r = gvalue_recorder.record_gvalue(gvalues)
fwdpy11.evolvets(rng, pop, params, 100, r)

assert len(gvalues) == pop.generation, "Callback failure"
