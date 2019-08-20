"""
Plots the mean genetic value for all time periods
present in the tree sequence.

Usage: python3 plot_genetic_values_from_tree_sequence.py filename

filename should contain a pickled population output by SlocusPopGSSmo.py
"""
import pickle
import lzma
import sys
import numpy as np
import seaborn as sns

with lzma.open(sys.argv[1], 'rb') as f:
    pop = pickle.load(f)

nodes = np.array(pop.tables.nodes, copy=False)

ancient_sample_metadata = np.array(pop.ancient_sample_metadata, copy=False)
alive_sample_metadata = np.array(pop.diploid_metadata, copy=False)
metadata = np.hstack((ancient_sample_metadata, alive_sample_metadata))

metadata_times = nodes['time'][metadata['nodes'][:, 0]]

mean_trait_values = []
for u in np.unique(metadata_times):
    samples_at_time_u = np.where(metadata_times == u)[0]
    mean_trait_values.append(
        metadata['g'][samples_at_time_u].mean())

# plot the time vs trait value
xyplot = sns.lineplot(np.unique(metadata_times), mean_trait_values)
xyplot.set_xlabel("Time (generations")
xyplot.set_ylabel("Mean trait value")
xyplot.get_figure().savefig("mean_trait_value_from_tree_sequences.png",
                            dpi=300)
