(migtest)=

# Symmetric migration between two demes

:::{figure-md}
<img src="../images/migtest.png">

Example output from this script, which simulates two demes
with `N=1,000` per deme.  For each replicate, {math}`F_{st}`
is calculated from a random sample of 10 diploids taken
from each deme at the end of the simulation. For the scaling
used here, the {math}`E[F_{st}]=\frac{1}{1+8Nm}`.  It took
about 17 minutes to get these results on a desktop machine
with an AMD Ryzen 3600 CPU.

:::

The script executes the following work flow:

* Determine a set of evenly-spaced migration rates from
  `m = 1e-5` up to `4Nm = 2`.
* For each migration rate, execute some number of replicate
  simulations in parallel using `concurrent.futures`.
  There are two demes, each with `N = 1,000`.
  The scaled recombination rate is {math}`\rho = 100` and migration
  is symmetrical between demes.
* For each replicate, add neutral mutations at scaled
  rate {math}`\theta = 500` and create an instance
  of {class}`fwdpy11.DataMatrix` based on ten
  randomly-chosen diploids from each deme.
* Calculate {math}`F_{st}` for each replicate
* Generate a plot like that shown above and save it to a file.

The `runsim` function handles the mechanics of the simulation itself.
The rest of the code is interface, execution, and plotting.

```{literalinclude} ../../examples/discrete_demography/migtest.py

```


