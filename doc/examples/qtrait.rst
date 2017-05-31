.. _qtraits1:

Simulating quantitative traits, I.
==========================================

Background reading:

* :ref:`definitions`

Further reading:

* :ref:`heritability`

fwdpy11 allows the simulation of quantitative traits in a fairly general way.  These simulations differ from standard
population genetic simulations in that a "trait value" is calculated and that a mapping from trait value to fitness is
required. Optionally, one may want to add "noise" to trait values, perhaps representing non-genetic contributions to
phenotype.

This page covers simulations of quantitative traits in a single genomic region (*e.g.* using
:class:`fwdpy11.fwdpy11_types.SlocusPop`).

Let's call our trait value :math:`P`, which is composed of a genetic component, :math:`G` and a random component
:math:`E`.

The evolve function
-----------------------------

The function used to evolve quantitative traits in :func:`fwdpy11.wright_fisher_qtrait.evolve_regions_sampler_fitness`.

Trait values
-----------------------------

The following trait value functions are implemented:

* :class:`fwdpy11.trait_values.SlocusAdditiveTrait`
* :class:`fwdpy11.trait_values.SlocusMultTrait`
* :class:`fwdpy11.trait_values.SlocusGBRTrait`

The above functions calculate the :math:`G` component of :math:`P`.  Custom trait value functions can be written in C++
in a similar manner as for custom fitness functions for standard "popgen" simulations.

Adding noise to trait values
----------------------------------------------------------

Let's look at the :math:`E` bit of :math:`P`.  The simplest scenario is :math:`P=G+E`, where :math:`E\sim N(\mu,\sigma)`.
Such a model is implemented in :class:`fwdpy11.wright_fisher_qtrait.GaussianNoise`:

.. literalinclude:: ../../fwdpy11/wright_fisher_qtrait.py
    :language: python
    :lines: 85-103

The implementation is straightforward:

* The class is constructed with a random number generator, the mean, and the standard deviation.
* The call operator gets passed three things: the offspring's :math:`G` value, and two objects of type
  :class:`fwdpy11.fwdpy11_types.SingleLocusDiploid` representing the two parents.  The function returns a Gaussian
  deviate with the appropriate mean and standard deviation.

For this example, the offspring genetic value and the parents are not needed. Below, we'll implement a simple model
where the offspring inherits the mean :math:`E` terms from its parents.

Mapping trait values to fitness
----------------------------------------------------------

In order to calulate fitness, we need a Python callable relating trait values to fitness.  In this section, I show the
implementation of two classes included in fwdpy11.

Stabilizing selection with constant parameters
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The first class models Gaussian stabilizing selection ("GSS") with respect to a constant optimum trait value and constant
intensity of selection. Under GSS models, :math:`w=e^{-\frac{(P-O)^2}{2VS}}`, where :math:`P` is the trait value,
:math:`O` is the optimum, and :math:`VS` determines the intensity of selection against deviations from the optimum
(*e.g* selection against variance in trait values). The class is :class:`fwdpy11.wright_fisher_qtrait.GSS`, and its
implementation is:

.. literalinclude:: ../../fwdpy11/wright_fisher_qtrait.py
    :language: python
    :lines: 12-38

This class is rather simple:

* The constructor takes the optimum and :math:`VS` as parameters.
* The call function takes a diploid's trait value and uses the class data to return :math:`w=e^{-\frac{(P-O)^2}{2VS}}`.

Stabilizing selection with a moving optimum
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The following code block shows the implementation of :class:`fwdpy11.wright_fisher_qtrait.GSSmo`, which models GSS with
changing optimum conditions:

.. literalinclude:: ../../fwdpy11/wright_fisher_qtrait.py
    :language: python
    :lines: 40-83

Again, this class is relatively simple:

* The constructor takes a list of tuples specifying the generations at which :math:`O` and/or :math:`VS` change. 
* The call operator is the same as in the previous example
* The new thing is the "update" function, which gets passed in the population's current generation.  We use the
  generation to check if the conditions change. If they do change, we update the model parameters and pop an element off
  the front of the list.

.. note::
    The "update" function is **optional**, and you only need to provide one if your trait value to fitness mapping
    function is "stateful" in some way, like :class:`fwdpy11.wright_fisher_qtrait.GSSmo`.

An example simulation
-----------------------------

The following code block represents the following model:

* There are no neutral mutations
* Mutations affecting trait values are additive and Gaussian distributed.
* Trait values are under Gaussian stabilizing selection with :math:`VS=1`.
* The optimum changes to a new value every :math:`0.1N` generations according to a geometric distribution.  The new optimum value is :math:`N(0,1)`.
* The :math:`E` value for offspring is the parental mean :math:`E` plus a Gaussian deviate.  Thus, there are shared
  environmental effects between parent and offspring that have a mean of zero.

.. testcode::

    import fwdpy11 as fp11
    import fwdpy11.trait_values as fp11tv
    import fwdpy11.wright_fisher_qtrait as fp11qt
    import fwdpy11.gsl_random as gsl
    import numpy as np

    class SharedE:
        def __init__(self,rng,sd):
            self.rng=rng
            self.sd=sd
            if(sd<0):
                raise ValueError("sd > 0 required")
        def __call__(self,p1,p2):
            mp = (p1.e+p2.e)/2.
            return mp + gsl.gsl_ran_gaussian_ziggurat(self.rng,self.sd)

    N = 1000
    pop = fp11.SlocusPop(N)


    rng = fp11.GSLrng(42)

    #Set up data for GSSmo
    gss_params=[(0,0.0,1.0)]
    #Starting at N gens into the sim,
    #the environment changes according to geometric
    #distribution with mean 0.1*N
    timepoint=N
    while timepoint <= 2*N:
        nt = gsl.gsl_ran_geometric(rng,1.0/float(0.1*N))
        timepoint += nt
        gss_params.append((timepoint,gsl.gsl_ran_gaussian_ziggurat(rng,1.0),1.0))

    t2f = fp11qt.GSSmo(gss_params)

    p = {'nregions':[],
    'sregions':[fp11.GaussianS(0,1,1,0.25)],
    'recregions':[fp11.Region(0,1,1)],
    'rates':(0.0,2e-3,1e-3),
    'demography':np.array([N]*N,dtype=np.uint32),
    'gvalue':fp11tv.SlocusAdditiveTrait(2.0),
    'trait_to_fitness':t2f,
    'noise':SharedE(rng,0.1)
    }

    params = fp11.model_params.SlocusParamsQ(**p)

    fp11qt.evolve(rng,pop,params)

