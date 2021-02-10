
Variation in mutation rates
=================================================================

Variation in mutation rates to neutral variants
------------------------------------------------------------------------------------------------

When simulating tree sequences, fwdpy11 currently does not allow variation in the neutral mutation rate.
Allowing this will come in future releases.

Variation in mutation rates to selected (non-neutral) variants
------------------------------------------------------------------------------------------------

Variation in the rate to non-neutral mutations is handled by classes derived from the ABC
:class:`fwdpy11.Sregion`.  Given an overall mutation rate to non-neutral variants, instances
of "sregions" are used to set up a multinomial distribution for generating new mutations.  

The following sets up a model where mutations have a constant effect size (:math:`s=-0.01`),
dominance :math:`h=0.25`, and occur uniformly on the interval :math:`[0, 1)`:

.. ipython:: python

    import fwdpy11

    sregions = [fwdpy11.ConstantS(beg=0, end=1, weight=1, s=-0.01, h=0.25)]
The previous example uses argument names for clarity, and the following is equivalent:

.. ipython:: python

    sregions = [fwdpy11.ConstantS(0, 1, 1, -0.01, 0.25)]
    print(sregions[0])


Note that the constructor parameters for these classes often have default values--see the specific class documentation 
for details.

In some scenarios, it is useful to think about the distribution of effect sizes as scaled with respect to the population
size.  For example, selection coefficients may be exponentially-distributed with a mean of :math:`2Ns`.  To do this in
fwdpy11:

.. ipython:: python

    # ALPHA = 2Ns
    MEAN_ALPHA = -10
    N = 1000
    sregions = [fwdpy11.ExpS(0, 1, 1, MEAN_ALPHA, scaling=2 * N)]

    print(sregions[0])


Region weighting
----------------------------------------

When multiple "sregion" objects are used, the default behavior is to multiply the input `weight` by `end-beg`:

.. ipython:: python

   sregions = [fwdpy11.ExpS(0, 1, 1, -0.2), fwdpy11.ConstantS(1, 3, 1, -0.1)]
   print(sregions)


Here, the input `weight` is interpreted to mean the weight "per site" is constant. In this example, twice as many mutations will have positions in :math:`[1, 3)` as from :math:`[0, 1)`.

To change the default behavior, one can prevent the coupling between input `weight` and region length:

.. ipython:: python

   sregions = [
       fwdpy11.ExpS(0, 1, 1, -0.2, coupled=False),
       fwdpy11.ConstantS(1, 3, 1, -0.1, coupled=False),
   ]
   print(sregions)


The absolute values of the `weight` parameters themselves is irrelevant.  The only thing that matters is the *relative* values from region to region.  Simulations based on the above examples would give the same results if the `weight` were 42 or 73.06.  Therefore, we can recreate our first example with code like the following:

.. ipython:: python

   sregions = [
       fwdpy11.ExpS(0, 1, 56, -0.2, coupled=False),
       fwdpy11.ConstantS(1, 3, 112, -0.1, coupled=False),
   ]
   print(sregions)


In the above example, twice as many mutations occur in the second region because the weights have relative values of 2:1.

.. note::

    Different regions are allowed to overlap, allowing the simulation of concepts like "coding regions"
    where the DFE are a weighted mixture from multiple distributions, etc.

