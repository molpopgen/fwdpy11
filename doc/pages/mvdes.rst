.. _mvdes:

.. ipython:: python
   :suppress:

   import fwdpy11
   import numpy as np

Different effect sizes of mutations in different demes
======================================================================

.. versionadded:: 0.7.0

This section describes how to allow mutations to have different effect sizes in different demes,
building on the material from :ref:`softselection`.

To allow for correlated effect sizes across demes, we used :class:`fwdpy11.mvDES` to model multivariate
distributions of effect sizes.  This class inherits from :class:`fwdpy11.Sregion` (see :ref:`mutationregions`).
In general, there is no standard *general* method for drawing random deviates from arbitrary multivariate
distributions.  The approach that we take is to use a multivariate Gaussian distribution as the underlying kernel.

First, we will describe the cases where the multivariate Gaussian kernel leads to output effect size distributions
that have straightforward interpretations.  Then, we will move onto the more general case allowing you to construct
multivariate distributions from current :class:`fwdpy11.Sregion` types.

The multivariate Gaussian
------------------------------------------------------------------

We may model Gaussian effect sizes using the existing :class:`fwdpy11.MultivariateGaussianEffects`
in conjunction with :class:`fwdpy11.mvDES`.  Using :class:`fwdpy11.MultivariateGaussianEffects` on its
own is used to model pleiotropic effects on a trait.  Here, we are re-using this type to model correlated
effect sizes across demes.

At this time, it is probably best to look at an example. The following code models Gaussian stabilizing
selection on a quantitative trait.  The effects sizes within each deme are themselves given by Gaussian
distributions and there is no correlation in the effect size in the two demes.

.. ipython:: python

    pdict = {'nregions': [],
            'recregions': [],
            'sregions': [fwdpy11.mvDES(fwdpy11.MultivariateGaussianEffects(0, 1, 1, np.identity(2)),
                         np.zeros(2))],
            'rates': (0, 1e-3, None),
            'demography': fwdpy11.DiscreteDemography(
                migmatrix=np.array([0.9, 0.1, 0.1, 0.9]).reshape((2, 2))
            ),
            'simlen': 100,
            'gvalue': fwdpy11.Additive(ndemes=2, scaling=2,
                                       gv2w=fwdpy11.GSS(opt=0, VS=10)),
           }

Most of the above is standard.  Let's dissect the new bits:

* An instance of :class:`fwdpy11.mvDES` is our only region with selected mutations.
* This instance holds an instance of :class:`fwdpy11.MultivariateGaussianEffects`
  that puts mutations on the interval :math:`[0, 1)` with weight 1 and an identity
  matrix specifies the correlation in effect sizes between demes 0 and 1.  The
  identity matrix has the value zero for all off-diagonal elements, meaning
  no covariance in effect sizes across demes.
* The final constructor argument specifies the mean of each marginal Gaussian
  distribution. The means are both zero.

Let's evolve the model now:

.. ipython:: python

    params = fwdpy11.ModelParams(**pdict)
    pop = fwdpy11.DiploidPopulation([100, 100], 1.0)
    rng = fwdpy11.GSLrng(1010)
    fwdpy11.evolvets(rng, pop, params, 10)

Let's extract the effect sizes from each deme:

.. ipython:: python

    for i in pop.tables.mutations:
        print(pop.mutations[i.key].esizes)

Let's look at another example where effect sizes covary negatively across demes and raise the mutation rate a bit:

.. ipython:: python

    vcv = np.array([1., -0.99, -0.99, 1.]).reshape((2,2))
    params.sregions = [fwdpy11.mvDES(fwdpy11.MultivariateGaussianEffects(0, 1, 1, vcv), np.zeros(2))]
    params.rates = (0, 5e-3, None)
    pop = fwdpy11.DiploidPopulation([100, 100], 1.0)
    fwdpy11.evolvets(rng, pop, params, 10)
    for i in pop.tables.mutations:
        print(pop.mutations[i.key].esizes)

Now we see that the effect sizes often differ in sign between the two demes.

The multivariate lognormal
------------------------------------------------------------------

If :math:`X` is a multivariate Gaussian distribution, :math:`N(\mathbf{\mu}, \mathbf{\sum})`, where :math:`\mathbf{\mu}` is a vector of mean values and 
:math:`\mathbf{\sum}` is the covariance matrix, then :math:`Y = e^X` is a
multivariate lognormal random variable with mean :math:`E[Y]_i = e^{\mu_i + \frac{1}{2}\sum_{ii}}` and covariance matrix :math:`Var[Y]_{i,j} = e^{\mu_i + \mu_j + \frac{1}{2}(\sum_{ii} + \sum_{jj})}(e^{\sum_{ij}}-1)`.

To specify a multivariate lognormal distribution of effect sizes, we use
the static class method :func:`fwdpy11.LogNormalS.mv`.  The following code
constructs a distribution of effect sizes such that `-2Ns` (where `N` is the 
size of a single deme) is a multivariate lognormal with means zero and an
identity matrix as a covariance matrix used to specify the multivate 
Gaussian kernel.

.. ipython:: python

    mvdes = fwdpy11.mvDES(fwdpy11.LogNormalS.mv(0, 1, 1, scaling=-200),
                np.zeros(2), np.identity(2))

.. note::

    The lognormal distribution returns deviates :math:`> 0`.
    To model deleterious mutations/effect sizes < 0, use the
    `scaling` parameter with a negative value like we just did!

Let's put it in a simulation and run it:

.. ipython:: python

    pdict = {'nregions': [],
            'recregions': [],
            'sregions': [mvdes],
            'rates': (0, 1e-3, None),
            'demography': fwdpy11.DiscreteDemography(
                migmatrix=np.array([0.9, 0.1, 0.1, 0.9]).reshape((2, 2))
            ),
            'simlen': 100,
            'gvalue': fwdpy11.Multiplicative(ndemes=2, scaling=2)
           }
    params = fwdpy11.ModelParams(**pdict)
    pop = fwdpy11.DiploidPopulation([100, 100], 1.0)
    fwdpy11.evolvets(rng, pop, params, 10)
    for i in pop.tables.mutations:
        print(pop.mutations[i.key].esizes)

"Custom" multivariate distributions
------------------------------------------------------------------

Recipes
------------------------------------------------------------------
