(recapitation)=

# Finishing a simulation with a tree sequence from msprime

:::{note}

See {ref}`here <msprime-subtleties>` for additional important information about interacting iwth `msprime`.

:::

This example is complementary to {ref}`precapitation`.  Rather than starting
with a tree sequence from `msprime`, we instead finish a simulation by "coalescing
back" the first generation of the simulation using `msprime`.  
{cite}`Haller2019-mn` refer to this procedure as "recapitation" of a tree sequence. In order for recapitation
to work correctly, we must pass `preserve_first_generation=True` to {func}`fwdpy11.evolvets`.
If we keep this option at its default value of `False` then we risk getting biased results.

The interface is the same as for the other example:

```{code-block} bash

time PYTHONPATH=../.. python3 recapitation.py -N 1000 --rho 1000 --theta 1000 --nreps 1000 --simlen 1000 --seed 42 --nsam 10 --model dtwf

```

The output format is also identical:

```
Means from fwdpy11:
S2N            8271.088000
Pi2N            998.004104
Sn             2823.524000
Pin             998.227800
watterson_n     998.075534
dtype: float64
Means from msprime:
S2N            8281.812000
Pi2N            999.956051
Sn             2827.652000
Pin             999.524778
watterson_n     999.534723
dtype: float64

real        8m25.782s
user        98m47.845s
sys 0m7.061s
```

:::{note}

The faster execution time of this example compared to {ref}`precapitation` is
*not* due to a performance difference due to starting with `msprime` vs finishing
with `msprime`.  Those two options are a wash in terms of performance.
The difference is that all of the statistics in this example are done using the
tree sequence statistics of `tskit` described in {cite}`Ralph2020-cg`. The `fwdpy11`
method {func}`fwdpy11.TableCollection.fs` is not (yet) implemented as a tree statistic
and calls to that function account for the bulk of the run time difference.

:::

Let's take a look at where the performance different is lurking.  If we execute the following
`iPython` block, we can profile the amount of time spent to setup and run the simulation
and calculate the summary statistics.

```{code-block} ipython

In [1]: import precapitate
   ...: import recapitate
   ...: %load_ext line_profiler
   ...: %lprun -f recapitate.run_sim recapitate.run_sim(1000, 1000, 1000, 1000, "dtwf", 100, 101, 666)
   ...: %lprun -f precapitate.run_sim precapitate.run_sim(1000, 1000, 1000, 1000, "dtwf", 100, 101, 666)
   ...:

```

The results for recapitation are:

```
Timer unit: 1e-06 s

Total time: 1.9008 s
File: /home/kevin/src/fwdpy11/examples/tskit/recapitate.py
Function: run_sim at line 86

Line #      Hits         Time  Per Hit   % Time  Line Contents
==============================================================
    86                                           def run_sim(N, rho, theta, simlen, model, nsam, fwdpy11_seed, msprime_seed):
    87         1         17.0     17.0      0.0      np.random.seed(msprime_seed)
    88         1         68.0     68.0      0.0      pop = fwdpy11.DiploidPopulation(N, 1.0)
    89
    90                                               pdict = {
    91         1          0.0      0.0      0.0          "nregions": [],
    92         1          1.0      1.0      0.0          "sregions": [],
    93                                                   "recregions": [
    94         1         16.0     16.0      0.0              fwdpy11.PoissonInterval(0, pop.tables.genome_length, rho / pop.N / 4)
    95                                                   ],
    96         1          0.0      0.0      0.0          "rates": (0, 0, None),
    97         1          5.0      5.0      0.0          "gvalue": fwdpy11.Multiplicative(1),
    98         1          5.0      5.0      0.0          "demography": fwdpy11.DiscreteDemography(),
    99         1          1.0      1.0      0.0          "simlen": simlen,
   100                                               }
   101         1         66.0     66.0      0.0      params = fwdpy11.ModelParams(**pdict)
   102         1         13.0     13.0      0.0      rng = fwdpy11.GSLrng(fwdpy11_seed)
   103         1     557747.0 557747.0     29.3      fwdpy11.evolvets(rng, pop, params, 100, preserve_first_generation=True)
   104         1      33930.0  33930.0      1.8      ts = pop.dump_tables_to_tskit()
   105         1    1261655.0 1261655.0     66.4      completed_ts = run_msprime(N, rho, 0.0, model, msprime_seed, ts)
   106         1        349.0    349.0      0.0      samples = pop.alive_nodes
   107         1          2.0      2.0      0.0      completed_ts = msprime.mutate(
   108         1      28703.0  28703.0      1.5          completed_ts, rate=theta / N / 4, random_seed=msprime_seed
   109                                               )
   110         1       4866.0   4866.0      0.3      S2N = completed_ts.segregating_sites([samples])[0]
   111         1       4440.0   4440.0      0.2      pi2N = completed_ts.diversity([samples])[0]
   112         1        100.0    100.0      0.0      rsample = np.random.choice(samples, nsam, replace=False)
   113         1       4415.0   4415.0      0.2      Sn = completed_ts.segregating_sites([rsample])[0]
   114         1       4398.0   4398.0      0.2      pin = completed_ts.diversity([rsample])[0]
   115
   116         1          6.0      6.0      0.0      return SimOutput(S2N, pi2N, Sn, pin)
```

The results for precapitation are:

```
Timer unit: 1e-06 s

Total time: 2.30427 s
File: /home/kevin/src/fwdpy11/examples/tskit/precapitate.py
Function: run_sim at line 96

Line #      Hits         Time  Per Hit   % Time  Line Contents
==============================================================
    96                                           def run_sim(N, rho, theta, simlen, model, nsam, fwdpy11_seed, msprime_seed):
    97         1         12.0     12.0      0.0      np.random.seed(msprime_seed)
    98         1    1187728.0 1187728.0     51.5      ts = run_msprime(N, rho, 0.0, model, msprime_seed)
    99         1       4168.0   4168.0      0.2      pop = fwdpy11.DiploidPopulation.create_from_tskit(ts)
   100         1          7.0      7.0      0.0      del ts  # No longer needed
   101
   102                                               pdict = {
   103         1          1.0      1.0      0.0          "nregions": [],
   104         1          0.0      0.0      0.0          "sregions": [],
   105                                                   "recregions": [
   106         1         20.0     20.0      0.0              fwdpy11.PoissonInterval(0, pop.tables.genome_length, rho / pop.N / 4)
   107                                                   ],
   108         1          1.0      1.0      0.0          "rates": (0, 0, None),
   109         1          6.0      6.0      0.0          "gvalue": fwdpy11.Multiplicative(1),
   110         1          7.0      7.0      0.0          "demography": fwdpy11.DiscreteDemography(),
   111         1          2.0      2.0      0.0          "simlen": simlen,
   112                                               }
   113         1         88.0     88.0      0.0      params = fwdpy11.ModelParams(**pdict)
   114         1          7.0      7.0      0.0      rng = fwdpy11.GSLrng(fwdpy11_seed)
   115         1     604449.0 604449.0     26.2      fwdpy11.evolvets(rng, pop, params, 100)
   116         1       4912.0   4912.0      0.2      fwdpy11.infinite_sites(rng, pop, theta / pop.N / 4)
   117         1        351.0    351.0      0.0      an = pop.alive_nodes
   118         1     251575.0 251575.0     10.9      S2N, pi2N = process_sim(pop, an, check_fixations=True)
   119         1         74.0     74.0      0.0      rn = np.random.choice(an, size=nsam, replace=False)
   120         1     250859.0 250859.0     10.9      Sn, pin = process_sim(pop, rn, check_fixations=False)
   121         1          4.0      4.0      0.0      return SimOutput(S2N, pi2N, Sn, pin)
```

It takes some squinting at first if you aren't used to looking at these line profiles.  Overall, the timing differences
between the two approaches are negligible for all of the steps except calculating the summary statistics. The
recapitation example is about 6 times faster due to the use of tree sequence statistics.

The main point is that precapitation and recapitation seem to be equally efficient.  The summary statistic calculation
is a bit of a distraction here, but it is perhaps useful to see that larger simulations make our (currently)
inefficient `fs` method less painful:

```{code-block} ipython

In [1]: import precapitate
   ...: import recapitate
   ...: %load_ext line_profiler
   ...: %lprun -f recapitate.run_sim recapitate.run_sim(5000, 5000, 5000, 5000, "dtwf", 100, 101, 666)
   ...: %lprun -f precapitate.run_sim precapitate.run_sim(5000, 5000, 5000, 5000, "dtwf", 100, 101, 666)
   ...:

```

```
Timer unit: 1e-06 s

Total time: 55.0555 s
File: /home/kevin/src/fwdpy11/examples/tskit/recapitate.py
Function: run_sim at line 86

Line #      Hits         Time  Per Hit   % Time  Line Contents
==============================================================
    86                                           def run_sim(N, rho, theta, simlen, model, nsam, fwdpy11_seed, msprime_seed):
    87         1         16.0     16.0      0.0      np.random.seed(msprime_seed)
    88         1        328.0    328.0      0.0      pop = fwdpy11.DiploidPopulation(N, 1.0)
    89
    90                                               pdict = {
    91         1          1.0      1.0      0.0          "nregions": [],
    92         1          1.0      1.0      0.0          "sregions": [],
    93                                                   "recregions": [
    94         1         16.0     16.0      0.0              fwdpy11.PoissonInterval(0, pop.tables.genome_length, rho / pop.N / 4)
    95                                                   ],
    96         1          0.0      0.0      0.0          "rates": (0, 0, None),
    97         1          5.0      5.0      0.0          "gvalue": fwdpy11.Multiplicative(1),
    98         1          6.0      6.0      0.0          "demography": fwdpy11.DiscreteDemography(),
    99         1          2.0      2.0      0.0          "simlen": simlen,
   100                                               }
   101         1         69.0     69.0      0.0      params = fwdpy11.ModelParams(**pdict)
   102         1         11.0     11.0      0.0      rng = fwdpy11.GSLrng(fwdpy11_seed)
   103         1   16544416.0 16544416.0     30.1      fwdpy11.evolvets(rng, pop, params, 100, preserve_first_generation=True)
   104         1     183217.0 183217.0      0.3      ts = pop.dump_tables_to_tskit()
   105         1   37901927.0 37901927.0     68.8      completed_ts = run_msprime(N, rho, 0.0, model, msprime_seed, ts)
   106         1        395.0    395.0      0.0      samples = pop.alive_nodes
   107         1          2.0      2.0      0.0      completed_ts = msprime.mutate(
   108         1     230681.0 230681.0      0.4          completed_ts, rate=theta / N / 4, random_seed=msprime_seed
   109                                               )
   110         1      46537.0  46537.0      0.1      S2N = completed_ts.segregating_sites([samples])[0]
   111         1      50670.0  50670.0      0.1      pi2N = completed_ts.diversity([samples])[0]
   112         1        229.0    229.0      0.0      rsample = np.random.choice(samples, nsam, replace=False)
   113         1      46352.0  46352.0      0.1      Sn = completed_ts.segregating_sites([rsample])[0]
   114         1      50568.0  50568.0      0.1      pin = completed_ts.diversity([rsample])[0]
   115
   116         1          8.0      8.0      0.0      return SimOutput(S2N, pi2N, Sn, pin)
Timer unit: 1e-06 s

Total time: 55.9958 s
File: /home/kevin/src/fwdpy11/examples/tskit/precapitate.py
Function: run_sim at line 96

Line #      Hits         Time  Per Hit   % Time  Line Contents
==============================================================
    96                                           def run_sim(N, rho, theta, simlen, model, nsam, fwdpy11_seed, msprime_seed):
    97         1         12.0     12.0      0.0      np.random.seed(msprime_seed)
    98         1   35346846.0 35346846.0     63.1      ts = run_msprime(N, rho, 0.0, model, msprime_seed)
    99         1      28786.0  28786.0      0.1      pop = fwdpy11.DiploidPopulation.create_from_tskit(ts)
   100         1         13.0     13.0      0.0      del ts  # No longer needed
   101
   102                                               pdict = {
   103         1          1.0      1.0      0.0          "nregions": [],
   104         1          1.0      1.0      0.0          "sregions": [],
   105                                                   "recregions": [
   106         1         22.0     22.0      0.0              fwdpy11.PoissonInterval(0, pop.tables.genome_length, rho / pop.N / 4)
   107                                                   ],
   108         1          0.0      0.0      0.0          "rates": (0, 0, None),
   109         1          6.0      6.0      0.0          "gvalue": fwdpy11.Multiplicative(1),
   110         1          5.0      5.0      0.0          "demography": fwdpy11.DiscreteDemography(),
   111         1          1.0      1.0      0.0          "simlen": simlen,
   112                                               }
   113         1         92.0     92.0      0.0      params = fwdpy11.ModelParams(**pdict)
   114         1          7.0      7.0      0.0      rng = fwdpy11.GSLrng(fwdpy11_seed)
   115         1   17592887.0 17592887.0     31.4      fwdpy11.evolvets(rng, pop, params, 100)
   116         1      38644.0  38644.0      0.1      fwdpy11.infinite_sites(rng, pop, theta / pop.N / 4)
   117         1        402.0    402.0      0.0      an = pop.alive_nodes
   118         1    1490600.0 1490600.0      2.7      S2N, pi2N = process_sim(pop, an, check_fixations=True)
   119         1        167.0    167.0      0.0      rn = np.random.choice(an, size=nsam, replace=False)
   120         1    1497318.0 1497318.0      2.7      Sn, pin = process_sim(pop, rn, check_fixations=False)
   121         1          5.0      5.0      0.0      return SimOutput(S2N, pi2N, Sn, pin)
```

```{literalinclude} ../../examples/tskit/recapitate.py

```


