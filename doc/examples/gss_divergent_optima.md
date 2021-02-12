(gssdivergentoptima)=

# Adaptation to different optima with multivariate distribution of effect sizes

This example uses a multivariate Gaussian distribution of effect sizes to simulate correlated mutational effects
on a single trait in two demes. See {ref}`mvdes` for background.

Two demes with `N` diploids each evolve for `10N` generations around an optimum
trait value of zero.  Then, the optimum changes to `1` in deme `0` and `-1` in
deme `1`.

The covariance matrix of effect sizes defines a strong positive correlation between
effect sizes in each deme.  Thus, after the optimum shift, adaptive mutations in
deme `0` are deleterious in deme `1` and vice-versa.

The simulation runs for 200 generations past the optimum shift.  At that point,
the mean trait value and mean fitness for each deme are printed to the screen.

:::{note}

If you want to model mutations having the same effect size in all demes,
then setting up the `sregions` and `gvalue` as done in {ref}`gss_vignette` will
do the trick.  You'll just need to modify the demography in that example
to include multiple demes.

:::

```{literalinclude} ../../examples/gss_divergent_optima/gss_divergent_optima.py

```


