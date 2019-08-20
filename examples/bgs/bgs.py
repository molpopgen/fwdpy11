"""
Recreates line 4 of Table 1 from
Hudson, R. R., and N. L. Kaplan. 1995.
“Deleterious Background Selection with Recombination.”
Genetics 141 (4): 1605–17.

The physical layout of the genome follows their Figure 1.
The fact that the neutral region is non-recombining
is a major simplification for the analysis.  By definition,
a single tree describes the history of that region.  If that
were not the case, we would have to take the weighted
mean of summaries of trees, with tree lengths (relative
to total genome length) as the weights.
"""
import fwdpy11
import numpy as np
import sys
import concurrent.futures

GENOME_LENGTH = 1.0
R = 0.04
U = 0.08
N = 1600
# number of diploids to sample
NSAM = 10


def runsim(argtuple):
    seed = argtuple

    rng = fwdpy11.GSLrng(seed)

    pdict = {'gvalue': fwdpy11.Multiplicative(2.),
             'rates': (0., U/2., R),  # The U/2. is from their eqn. 2.
             'nregions': [],
             'sregions': [fwdpy11.ConstantS(0, 1./3., 1, -0.02, 1.),
                          fwdpy11.ConstantS(2./3., 1., 1, -0.02, 1.)],
             'recregions': [fwdpy11.Region(0, 1./3., 1),
                            fwdpy11.Region(2./3., 1., 1)],
             'demography': np.array([N]*20*N, dtype=np.uint32)
             }
    params = fwdpy11.ModelParams(**pdict)

    pop = fwdpy11.DiploidPopulation(N, GENOME_LENGTH)

    fwdpy11.evolvets(rng, pop, params, 100, suppress_table_indexing=True)

    rdips = np.random.choice(N, NSAM, replace=False)
    md = np.array(pop.diploid_metadata, copy=False)
    rdip_nodes = md['nodes'][rdips].flatten()
    nodes = np.array(pop.tables.nodes, copy=False)
    # Only visit trees spanning the
    # mutation-free segment of the genome
    tv = fwdpy11.TreeIterator(pop.tables, rdip_nodes, begin=1./3., end=2./3.)
    plist = np.zeros(len(nodes), dtype=np.int8)
    sum_pairwise_tmrca = 0
    for t in tv:
        for i in range(len(rdip_nodes)-1):
            u = rdip_nodes[i]
            while u != fwdpy11.NULL_NODE:
                plist[u] = 1
                u = t.parent(u)
            for j in range(i+1, len(rdip_nodes)):
                u = rdip_nodes[j]
                while u != fwdpy11.NULL_NODE:
                    if plist[u] == 1:
                        sum_pairwise_tmrca += 2 * \
                            (pop.generation-nodes['time'][u])
                        u = fwdpy11.NULL_NODE
                    else:
                        u = t.parent(u)
            plist.fill(0)
    return 2*sum_pairwise_tmrca/(len(rdip_nodes)*(len(rdip_nodes)-1))


if __name__ == "__main__":
    seed = int(sys.argv[1])
    nreps = int(sys.argv[2])
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
        futures = {executor.submit(runsim, (i)) for i in seeds}
        for future in concurrent.futures.as_completed(futures):
            scaled_mean_pairwise_coalescence_time = future.result()/(4*N)
            mean_pairwise_coalescence_times[idx] = scaled_mean_pairwise_coalescence_time
            idx += 1
    mean_pairwise_coalescence_times
    print(mean_pairwise_coalescence_times.mean(),
          mean_pairwise_coalescence_times.std()/np.sqrt(nreps))
