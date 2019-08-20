"""
"""
import pickle
import lzma
import sys
import numpy as np
import seaborn as sns
import fwdpy11

with lzma.open(sys.argv[1], 'rb') as f:
    pop = pickle.load(f)

nodes = np.array(pop.tables.nodes, copy=False)

ancient_sample_metadata = np.array(pop.ancient_sample_metadata, copy=False)
alive_sample_metadata = np.array(pop.diploid_metadata, copy=False)
metadata = np.hstack((ancient_sample_metadata, alive_sample_metadata))

metadata_nodes = metadata['nodes'].flatten()
metadata_node_times = nodes['time'][metadata_nodes]
metadata_record_times = nodes['time'][metadata['nodes'][:, 0]]


genetic_trait_values_from_sim = []
genetic_values_from_ts = []
idx = 0
for u in np.unique(metadata_node_times):
    samples_at_time_u = metadata_nodes[np.where(metadata_node_times == u)]
    vi = fwdpy11.VariantIterator(
        pop.tables, samples_at_time_u)
    sum_esizes = np.zeros(len(samples_at_time_u))
    for variant in vi:
        g = variant.genotypes
        r = variant.records[0]
        mutant = np.where(g == 1)[0]
        sum_esizes[mutant] += pop.mutations[r.key].s
    ind = int(len(samples_at_time_u)/2)
    temp_gvalues = np.zeros(ind)
    temp_gvalues += sum_esizes[0::2]
    temp_gvalues += sum_esizes[1::2]
    genetic_values_from_ts.extend(temp_gvalues.tolist())
    genetic_trait_values_from_sim.extend(
        metadata['g'][np.where(metadata_record_times == u)[0]].tolist())

# plot the time vs trait value
xyplot = sns.scatterplot(genetic_trait_values_from_sim, genetic_values_from_ts)
xyplot.set_xlabel("From simulation")
xyplot.set_ylabel("From tree sequences")
xyplot.get_figure().savefig("compare_genetic_values.png",
                            dpi=300)
