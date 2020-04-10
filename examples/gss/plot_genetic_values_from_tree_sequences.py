"""
Plots the mean genetic value for all time periods
present in the tree sequence.

Usage: python3 plot_genetic_values_from_tree_sequence.py filename

filename should contain a pickled population output by SlocusPopGSSmo.py
"""
import lzma
import pickle
import sys

import numpy as np
import seaborn as sns

with lzma.open(sys.argv[1], "rb") as f:
    pop = pickle.load(f)

nodes = np.array(pop.tables.nodes, copy=False)

times = []
mean_trait_values = []
for t, n, m in pop.sample_timepoints():
    times.append(t)
    mean_trait_values.append(m["g"].mean())

# plot the time vs trait value
xyplot = sns.lineplot(times, mean_trait_values)
xyplot.set_xlabel("Time (generations")
xyplot.set_ylabel("Mean trait value")
xyplot.get_figure().savefig("mean_trait_value_from_tree_sequences.png", dpi=300)
