.. _gss:

Gaussian stabilizing selection on a single trait (with a sudden optimum shift)
==============================================================================================================================

This simulation script shows how to simulate stabilizing selection on a single trait.  The fitness model is Gaussian
stabilizing selection.  The example is quite detailed, with a command-line interface implemented using `argparse`.

The simulation can track individuals as nodes in tree sequences after the time of the optimum shift.  For example,
to preserve 50 randomly-chosen individuals every 10 generations (starting from the time of the shift):

.. code-block:: bash

    python3 DiploidPopulationGSSmo.py -N 5000 --mu 5e-3 --sigma 0.1 --preserve 10 --num_ind 50 --filename pop.lzma

The output file, `pop.lzma`, is a compressed file containing the pickled population and its corresponding tree sequence.

.. literalinclude:: ../../examples/gss/DiploidPopulationGSSmo.py

Processing the samples stored in tree sequences
----------------------------------------------------------------------------------

This script plots the mean trait value at every time point as a function of time:

.. literalinclude:: ../../examples/gss/plot_genetic_values_from_tree_sequences.py

This script goes over every time point, iterates over the variants, and re-calculates the genetic values.
The output is a plot of the genetic value stored in the metadata versus the recalculated value.  The reason we can do 
this is that the simulation assumes strict additivity.

.. literalinclude:: ../../examples/gss/iterate_variants_in_tree_sequences.py

