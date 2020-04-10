"""
"""
import lzma
import pickle
import sys

import numpy as np
import seaborn as sns

import fwdpy11

with lzma.open(sys.argv[1], "rb") as f:
    pop = pickle.load(f)

genetic_trait_values_from_sim = []
genetic_values_from_ts = []
idx = 0
for n, s, m in pop.sample_timepoints():
    vi = fwdpy11.VariantIterator(pop.tables, s)
    sum_esizes = np.zeros(len(s))
    for variant in vi:
        g = variant.genotypes
        r = variant.records[0]
        mutant = np.where(g == 1)[0]
        sum_esizes[mutant] += pop.mutations[r.key].s
    ind = int(len(s) / 2)
    temp_gvalues = np.zeros(ind)
    temp_gvalues += sum_esizes[0::2]
    temp_gvalues += sum_esizes[1::2]
    genetic_values_from_ts.extend(temp_gvalues.tolist())
    genetic_trait_values_from_sim.extend(m["g"].tolist())

# plot the time vs trait value
xyplot = sns.scatterplot(genetic_trait_values_from_sim, genetic_values_from_ts)
xyplot.set_xlabel("From simulation")
xyplot.set_ylabel("From tree sequences")
xyplot.get_figure().savefig("compare_genetic_values.png", dpi=300)
