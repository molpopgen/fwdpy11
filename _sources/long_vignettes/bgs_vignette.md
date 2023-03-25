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

(bgs_vignette)=

# Background selection

This vignette recreates results from {cite}`Hudson1995-id`.

The physical layout of the genome follows their Figure 1.
The fact that the neutral region is non-recombining is a major simplification for the analysis.
By definition, a single tree describes the history of that region.
If that were not the case, we would have to take the weighted mean of summaries of trees, with tree lengths (relative to total genome length) as the weights.

```{code-cell} python
import numpy as np

import fwdpy11

def runsim(N, R, U, nsam, seed):
    rng = fwdpy11.GSLrng(seed)

    pdict = {
        "gvalue": fwdpy11.Multiplicative(2.0),
        "rates": (0.0, U / 2.0, None),  # The U/2. is from their eqn. 2.
        "nregions": [],
        "sregions": [
            fwdpy11.ConstantS(0, 1.0 / 3.0, 1, -0.02, 1.0),
            fwdpy11.ConstantS(2.0 / 3.0, 1.0, 1, -0.02, 1.0),
        ],
        "recregions": [
            fwdpy11.PoissonInterval(0, 1.0 / 3.0, R/2.0),
            fwdpy11.PoissonInterval(2.0 / 3.0, 1.0, R/2.0),
        ],
        # Evolve a single deme of size N for 20*N 
        # generations
        "demography": fwdpy11.ForwardDemesGraph.tubes([N], 20),
        "simlen": 20 * N,
    }
    params = fwdpy11.ModelParams(**pdict)

    pop = fwdpy11.DiploidPopulation(N, 1.0)

    fwdpy11.evolvets(rng, pop, params, 100, suppress_table_indexing=True)

    rdips = np.random.choice(N, nsam, replace=False)
    md = np.array(pop.diploid_metadata, copy=False)
    rdip_nodes = md["nodes"][rdips].flatten()
    nodes = np.array(pop.tables.nodes, copy=False)
    # Only visit trees spanning the
    # mutation-free segment of the genome
    tv = fwdpy11.TreeIterator(pop.tables, rdip_nodes, begin=1.0 / 3.0, end=2.0 / 3.0)
    plist = np.zeros(len(nodes), dtype=np.int8)
    sum_pairwise_tmrca = 0
    for t in tv:
        for i in range(len(rdip_nodes) - 1):
            u = rdip_nodes[i]
            while u != fwdpy11.NULL_NODE:
                plist[u] = 1
                u = t.parent(u)
            for j in range(i + 1, len(rdip_nodes)):
                u = rdip_nodes[j]
                while u != fwdpy11.NULL_NODE:
                    if plist[u] == 1:
                        sum_pairwise_tmrca += 2 * (pop.generation - nodes["time"][u])
                        u = fwdpy11.NULL_NODE
                    else:
                        u = t.parent(u)
            plist.fill(0)
    return 2 * sum_pairwise_tmrca / (len(rdip_nodes) * (len(rdip_nodes) - 1))
```

```{code-cell} python
---
tags: ['remove-cell']
---
# Dummy cell that doesn't get written to docs.
# This makes sure that we have fast examples
# that work.
R = 0.04
U = 0.01
N = 100
# number of diploids to sample
NSAM = 10
seed = 54321
runsim(N, R, U, NSAM, seed)
```

The following function runs the parameters in line 4 of Table 1 from {cite}`Hudson1995-id`.

```{code-cell} python
def hk95_table1_line4(seed):
    R = 0.04
    U = 0.08
    N = 1600
    # number of diploids to sample
    NSAM = 10
    return runsim(N, R, U, NSAM, seed), N
```

The following code block is not executed.
If you download this page as a notebook and execute it, you should get a mean diversity very close to $0.368$.

```{code-block} python
import concurrent.futures

nreps = 100
seed = 42
np.random.seed(seed)

seeds = []
for i in range(nreps):
    candidate = np.random.randint(0, np.iinfo(np.uint32).max)
    while candidate in seeds:
        candidate = np.random.randint(0, np.iinfo(np.uint32).max)
    seeds.append(candidate)

mean_pairwise_coalescence_times = np.zeros(nreps)
mean_pairwise_coalescence_times.fill(np.nan)
idx = 0
with concurrent.futures.ProcessPoolExecutor() as executor:
    futures = {executor.submit(hk95_table1_line4, i) for i in seeds}
    for future in concurrent.futures.as_completed(futures):
        scaled_mean_pairwise_coalescence_time, N = future.result() 
        scaled_mean_pairwise_coalescence_time /= (4 * N)
        mean_pairwise_coalescence_times[idx] = scaled_mean_pairwise_coalescence_time
        idx += 1
mean_pairwise_coalescence_times
print(
    mean_pairwise_coalescence_times.mean(),
    mean_pairwise_coalescence_times.std() / np.sqrt(nreps),
)
```
