Gaussian stabilizing selection in a single contiguous genomic interval
==========================================================================================

The scripts here involve simulating a continuous phenotype (without pleiotropy) along a 
continuous genomic segment where mutations occur as a uniform process on the half-open
interval [0,1).  

The simulations are run using tree sequence recording.  The output is a pickled population and,
optionally, an `sqlite3` database of quantitative-genetic summaries of the population over time.

The scripts to process the data work by processing either the `sqlite3` database or the tree sequences.
For database processing, I use `R`'s excellent dplyr_ package, which writes the SQL code for you! See here_
for more info about that.

Scripts
-------------------------------------

The simulation script is `SlocusPopGSSmo.py`.  To see how to use it:

.. code-block:: bash

    python3 SlocusPopGSSmo.py --help

There are several scripts:

* `plotstats.R` is an executable `R` script that plots the quantitative-genetic summaries recorded during the
  simulation, assuming that option was used.  It uses dplyr_ to read the summaries from an `sqlite3` database.
* `plot_genetic_values_from_tree_sequences.py` goes over time series samples recorded into tree sequences during the
  simulation.  For each time point, the mean genetic value is obtained and plotted over time.
* `iterate_variants_in_tree_sequences.py` also iterates over every time point on the tree sequences.  For each time
  point, the script iterates over the variants for each genotype and recalculates the genetic value for each individual.
  The final result in a plot of the genetic values recorded in the individual metadata versus the reconstructed values.
  The plot shows that `x = y` because the metadata are stored in order of the birth time of individuals, and the genetic
  values are easily reconstructed under the additive trait model simulated.  In effect, this script is a form of an
  integration test showing that the variants are properly traversed at each time point.

An example workflow
++++++++++++++++++++++++++++++++++++++++

The `Makefile` provided shows an example of using all of the above scripts.  To run it:

.. code-block:: bash

    make

.. _dplyr: https://dplyr.tidyverse.org/
.. _here: https://db.rstudio.com/dplyr/
