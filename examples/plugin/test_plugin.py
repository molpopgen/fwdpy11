import demes
import fwdpy11
import gvalue_recorder

N = 1000

yaml = """
time_units: generations
demes:
 - name: pop
   epochs:
    - start_size: 1000
"""

graph = demes.loads(yaml)
demography = fwdpy11.discrete_demography.from_demes(graph, 1)

pdict = {
    "demography": demography,
    "simlen": 100,
    "nregions": [],
    "sregions": [fwdpy11.GaussianS(0, 1, 1, 0.25)],
    "recregions": [fwdpy11.PoissonInterval(0, 1, 5e-3)],
    "rates": (0.0, 5e-3, None),
    "gvalue": fwdpy11.Additive(
        2.0,
        fwdpy11.GaussianStabilizingSelection.single_trait(
            optima=[fwdpy11.Optimum(VS=1.0, optimum=1.0, when=0)]
        ),
    ),
    "prune_selected": False,
}

params = fwdpy11.ModelParams(**pdict)

pop = fwdpy11.DiploidPopulation(N, 1.0)

rng = fwdpy11.GSLrng(42)

gvalues = []
rg = gvalue_recorder.record_gvalue()


def r(pop, _):
    rg(pop, gvalues)


fwdpy11.evolvets(rng, pop, params, 100, r)

assert len(gvalues) == pop.generation, "Callback failure"
