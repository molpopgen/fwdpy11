.. _startingfrommsprime:

Starting a simulation with a tree sequence from msprime
================================================================================

The following command line uses :mod:`msprime` to simulate under the discrete
time Wright-Fisher model using the methods described in [Nelson2019]_.  Then,
`fwdpy11` simulates for another 1,000 generations.  From the end of each simulation,
we get the total number of segregating sites in the entire population and a sample
of 10 random genomes.  We then do the same thing for a large number of independent
replicates from `msprime`.

.. code-block:: bash

    time PYTHONPATH=../.. python3 init_with_ts.py -N 1000 --rho 1000 --theta 1000 --nreps 1000 --simlen 1000 --seed 42 --nsam 10 --model dtwf

The output from the command is shown below and the timings are from a 12 core Intel Xeon W-2135 
processor:

::

    Means from fwdpy11:
    S2N            8275.392000
    Pi2N            998.949229
    Sn             2825.591000
    Pin             999.423778
    watterson_n     998.806189
    dtype: float64
    Means from msprime:
    S2N            8281.812000
    Pi2N            999.956051
    Sn             2828.331000
    Pin             999.089667
    watterson_n     999.774740
    dtype: float64

    real	9m30.625s
    user	110m18.927s
    sys	0m13.800s

.. literalinclude:: ../../examples/tskit/init_with_ts.py

