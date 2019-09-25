import numpy as np


def _get_times(pop):
    amd = np.array(pop.ancient_sample_metadata, copy=False)
    nodes = np.array(pop.tables.nodes, copy=False)
    times = nodes['time'][amd['nodes'][:, 0]]
    utimes = np.unique(times)
    return times, utimes


def traverse_sample_timepoints(pop, include_alive):
    amd = np.array(pop.ancient_sample_metadata, copy=False)
    amd.flags.writeable = False
    times, utimes = _get_times(pop)
    for uti in utimes:
        mdidx = np.where(times == uti)[0]
        sample_nodes = amd['nodes'][mdidx].flatten()
        sample_nodes.flags.writeable = False
        mdslice = amd[mdidx][:]
        mdslice.flags.writeable = False
        yield uti, sample_nodes, mdslice

    if include_alive is True:
        md = np.array(pop.diploid_metadata, copy=False)
        nodes = md['nodes'].flatten()
        md.flags.writeable = False
        nodes.flags.writeable = False
        yield pop.generation, nodes, md
