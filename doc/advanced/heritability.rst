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
    import fwdpy11.genetic_values
    import fwdpy11.genetic_value_noise
    import fwdpy11.model_params
    import math

.. ipython:: python

    rng = fwdpy11.GSLrng(42)
    # Parameters dict with some
    # arbitrary stuff in there:
    # Set VS = 1-VE
    gv2w = fwdpy11.genetic_values.GSS(VS=1-VE,opt=0)
    # Set sigma of noise function to square root of VE
    noise = fwdpy11.genetic_value_noise.GaussianNoise(mean=0,sd=np.sqrt(VE))
    p = {'nregions': [fwdpy11.Region(0,1,1)],
         'sregions': [fwdpy11.ExpS(0, 1, 1, 0.25)],
         'recregions': [fwdpy11.Region(0, 1, 1)],
         'rates': (1e-3, 2e-3, 1e-3),
         'gvalue': fwdpy11.genetic_values.DiploidAdditive(2.0,gv2w,noise),
         'prune_selected': False,
         }
    params = fwdpy11.model_params.ModelParams(**p)

Note that `gv2w` and `noise` would normally be written into `p` along with everything else.  However, Sphinx suppresses
comments in multi-line code chunks, requiring us to write the commands out separately from the rest of the dict.
