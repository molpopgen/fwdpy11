---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

(incomplete_selective_sweeps)=

# Incomplete selective sweeps

A sweep from a new mutation simulated until it first
hits a frequency :math:`\geq 0.25`.

Key points:

* Note the flags passed to {func}`fwdpy11.conditional_models.selective_sweep`.
* Note the differences in parameter scaling between this software and `msprime`.
* This example uses a custom monitor object, `IncompleteSweep`.

```{code-cell} python
import fwdpy11
import numpy as np
import msprime
import fwdpy11.conditional_models
import fwdpy11.tskit_tools


class IncompleteSweep(object):
    def __call__(
        self, pop: fwdpy11.DiploidPopulation, index: int, key: tuple
    ) -> fwdpy11.conditional_models.SimulationStatus:
        if pop.mutations[index].key != key:
            # it is fixed or lost, neither of 
            # which we want
            return fwdpy11.conditional_models.SimulationStatus.Restart
        if pop.mcounts[index] == 0:
            return fwdpy11.conditional_models.SimulationStatus.Restart
        
        # Terminate the first time we see the 
        # variant get about a freq of 0.25
        if pop.mcounts[index] / 2 / pop.N >= 0.25:
            return fwdpy11.conditional_models.SimulationStatus.Success
        # make sure there's a valid return value
        return fwdpy11.conditional_models.SimulationStatus.Continue

L = 10000.0
ttl_rec_rate = 1e-5*L

def setup(prune_selected=False):
    # Dropping mutations requires existing
    # ancestry, which we can get either
    # from a burn-in or from msprime.
    initial_ts = msprime.sim_ancestry(
        samples=500,
        population_size=500,
        # In msprime, recombination rates are per "unit" ("base pair")
        recombination_rate=ttl_rec_rate/L,
        random_seed=43215,
        sequence_length=L,
    )

    # Build the pop from msprime output
    pop = fwdpy11.DiploidPopulation.create_from_tskit(initial_ts)

    # Set up basic model parameters
    pdict = {
        # Here, the rec rate is number of events per region, not per "base pair"!
        "recregions": [fwdpy11.PoissonInterval(0, int(L), ttl_rec_rate, discrete=True)],
        "gvalue": fwdpy11.Multiplicative(2.0),
        "rates": (0, 0, None),
        "prune_selected": False,
        "simlen": 200,
    }
    params = fwdpy11.ModelParams(**pdict)

    return pop, params

ALPHA = 1000.0
rng = fwdpy11.GSLrng(12345)
pop, params = setup()
mutation_data = fwdpy11.conditional_models.NewMutationParameters(
    frequency=fwdpy11.conditional_models.AlleleCount(1),
    data=fwdpy11.NewMutationData(effect_size=ALPHA / 2 / pop.N, dominance=1),
    position=fwdpy11.conditional_models.PositionRange(left=L/2-1, right=L/2+1),
)
output = fwdpy11.conditional_models.selective_sweep(
    rng, 
    pop,
    params,
    mutation_data,
    IncompleteSweep(),
    # NOTE: this is key flag!
    # The default is False, which will 
    # keep simulating until params.simlen,
    # at which point the mutation may be fixed.
    return_when_stopping_condition_met = True
)


print(output.pop.mutations[output.mutation_index])

ts = output.pop.dump_tables_to_tskit()
print("genotype table of resultant sweep")
print(ts.genotype_matrix())
```
