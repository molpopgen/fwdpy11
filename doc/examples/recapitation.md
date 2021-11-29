(recapitation)=

# Finishing a simulation with a tree sequence from msprime

:::{note}

See {ref}`here <msprime-subtleties>` for additional important information about interacting with `msprime`.

:::

This example is complementary to {ref}`precapitation`.
Rather than starting with a tree sequence from `msprime`, we instead finish a simulation by "coalescing back" the first generation of the simulation using `msprime`.
{cite}`Haller2019-mn` refer to this procedure as "recapitation" of a tree sequence.
In order for recapitation to work correctly, we must pass `preserve_first_generation=True` to {func}`fwdpy11.evolvets`.
If we keep this option at its default value of `False` then we risk getting biased results.

The interface is the same as for the other example:

```{code-block} bash

time PYTHONPATH=../.. python3 recapitate.py -N 1000 --rho 1000 --theta 1000 --nreps 1000 --simlen 1000 --seed 42 --nsam 10 --model dtwf

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

```

```{literalinclude} ../../examples/tskit/recapitate.py

```


