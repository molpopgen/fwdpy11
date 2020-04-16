.. _finishwithmsprime:

Finishing a simulation with a tree sequence from msprime
================================================================================

.. code-block:: bash

    time PYTHONPATH=../.. python3 finish_with_msprime.py -N 1000 --rho 1000 --theta 1000 --nreps 1000 --simlen 1000 --seed 42 --nsam 10 --model dtwf

::

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

    real	8m25.782s
    user	98m47.845s
    sys	0m7.061s

.. literalinclude:: ../../examples/tskit/finish_with_msprime.py
