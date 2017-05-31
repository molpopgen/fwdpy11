.. _heritability:

Quantitative trait simulations with heritability
======================================================================

Background:

* :ref:`qtraits1`

Under a fairly wide range of conditions, simulations of quantitative traits depend largely on a single parameter, :math:`VS`.  This parameter can be broken down into :math:`VS = \omega^2 + VE`, where :math:`VE` is the non-genetic component to variation in trait values.  The expected genetic variance under a wide range of assumptions is approximately :math:`E[VG]\approx 4\mu VS`.  Thus, if you want to simulate a trait with a specific heritability (:math:`H^2 = VG/(VG+VE)`), you need to account for :math:`VE` when parameterizing a simulation.  When objects in this pacakge refer to :math:`VS`, I'm using a shorthand notation for :math:`\omega^2+VE`, meaing that if you set :math:`VS=1` and then add Gaussian noise to your trait such that :math:`VE \sim N(0,\sigma_E^2)`, then your simulation has :math:`VS=1+VE`, which will affect the average :math:`VG` in your output.

Let's work through a specific example.  First, define some parameters:

.. ipython:: python

    mu=1e-3
    VS=1.0
    h2=0.25
    EVG=4.*mu*VS
    print("E[VG] = ",EVG)

What should our :math:`VE` be?  That's easy to get by rearranging the usual formula for heritability:

.. ipython:: python

    VE = EVG*(1.0 - h2)/h2
    print("VE should be ",VE)

Thus, when we parameterize objects for our simulations, we should only pass :math:`\omega^2` in for :math:`VS`, and we should use our :math:`VE` to parameterize our noise function:

.. ipython:: python
    :suppress:

    import fwdpy11
    import fwdpy11.wright_fisher_qtrait 
    import math

.. ipython:: python

    rng = fwdpy11.GSLrng(42)
    trait_to_fitness = fwdpy11.wright_fisher_qtrait.GSS(VS-VE,0)
    noise = fwdpy11.wright_fisher_qtrait.GaussianNoise(rng,math.sqrt(VE))

