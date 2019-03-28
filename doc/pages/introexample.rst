Example
======================================================================


.. ipython:: python

    import fwdpy11
    import numpy as np

    N = 1000
    GENOME_LEN = 1.0
    THETA, RHO = 1000.0, 1000.0
    MU = 1e-3
    
    gssmo = fwdpy11.GSSmo([(0,0,1), (10*N,1,1)]) 
    pdict = {'gvalue': fwdpy11.Additive(2.0, gssmo),
            'nregions': [],
            'sregions': [fwdpy11.GaussianS(0, 1, 1, 0.15, 1)],
            'recregions': [fwdpy11.Region(0,1,1)],
            'rates': (0.0, MU, RHO/(4*N)),
            'demography': np.array([N]*(10*N + 100), dtype=np.uint32),
            'prune_selected': False
            }
    params = fwdpy11.ModelParams(**pdict)

    pop = fwdpy11.DiploidPopulation(N, GENOME_LEN)

    rng = fwdpy11.GSLrng(42)

    class Recorder(object):
        def __init__(self):
            self.gbar = []
        def __call__(self, pop, ancient_sampler_recorder):
            if pop.generation >= 10*pop.N:
                md = np.array(pop.diploid_metadata, copy=False)
                self.gbar.append((pop.generation, md['g'].mean()))
                ancient_sampler_recorder.assign(np.arange(pop.N, dtype=np.int32))

    recorder = Recorder()
    fwdpy11.evolvets(rng, pop, params, 100, recorder)

    print(pop.generation)

    # Let's get the mean trait value, the genetic variance and fitness
    # for the current generation
    alive_metadata = np.array(pop.diploid_metadata, copy=False)
    print(alive_metadata.dtype)
    print(alive_metadata['g'].mean(), alive_metadata['g'].var(), alive_metadata['w'].mean())

    print(np.mean([i.g for i in pop.diploid_metadata]))
    print(np.var([i.g for i in pop.diploid_metadata]))
    print(np.mean([i.w for i in pop.diploid_metadata]))


Plot the mean genetic value over time:

.. ipython:: python

    ancient_md = np.array(pop.ancient_sample_metadata, copy = False)
    node_table = np.array(pop.tables.nodes, copy=False)
    print(node_table.dtype)
    ancient_md_times = node_table['time'][ancient_md['nodes'][:,0]]

.. ipython:: python

    from matplotlib import rc
    rc('font',**{'size':18})
    rc('text', usetex=True)
    import matplotlib.pyplot as plt
    mean_genetic_values = []
    for t in np.unique(ancient_md_times):
        samples_at_t = np.where(ancient_md_times == t)[0]
        mean_genetic_values.append(ancient_md['g'][samples_at_t].mean())

    plt.plot(np.unique(ancient_md_times), mean_genetic_values);
    plt.ylabel("Mean trait value");
    plt.title("Adaptive walk to new optimum");
    plt.xlabel("Generation");
    @savefig mean_genetic_values_over_time.png width=6in
    plt.tight_layout();

Sanity check our calculations:

.. TODO::

    comment on np.concatenate to merge ancient + alive metadata

.. ipython:: python

    assert all([i==j[1] for i,j in zip(mean_genetic_values,recorder.gbar[:-1])]) is True
    assert recorder.gbar[-1][1] == alive_metadata['g'].mean()


.. ipython:: python

    all_md = np.concatenate((ancient_md, alive_metadata))
    print(all_md.flags)
    mdtimes = node_table['time'][all_md['nodes'][:,0]]
    ssh_over_time = []
    nmuts = fwdpy11.infinite_sites(rng, pop, THETA/(4*N))
    np.random.seed(54321)
    for t in np.unique(mdtimes):
        samples_at_t = np.where(mdtimes == t)[0]
        rsamples = np.random.choice(samples_at_t, 25, replace=False)
        rsamples_nodes = all_md['nodes'][rsamples,:].flatten()
        vi = fwdpy11.VariantIterator(pop.tables, pop.mutations, rsamples_nodes)
        ssh = 0.0
        for v in vi:
            g = v.genotypes
            r = v.record
            if pop.mutations[r.key].neutral is True:
                daf = float(g.sum())
                het = 2*daf*(len(g)-daf)/float(len(g)*(len(g)-1))
                ssh += het
        ssh_over_time.append(ssh)


    plt.plot(np.unique(mdtimes), ssh_over_time);
    plt.ylabel(r'$\pi$');
    @savefig pi_over_time.png width=6in
    plt.xlabel("Generation");
