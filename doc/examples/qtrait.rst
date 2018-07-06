.. _qtraits1:

Simulating quantitative traits, I
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
:class:`fwdpy11.SlocusPop`).

Let's call our trait value :math:`P`, which is composed of a genetic component, :math:`G` and a random component
:math:`E`.

The evolve function
-----------------------------

The function used to evolve quantitative traits in :func:`fwdpy11.wright_fisher`.

Trait values
-----------------------------

The following trait value functions are implemented:

* :class:`fwdpy11.genetic_values.SlocusAdditive`
* :class:`fwdpy11.genetic_values.SlocusMult`
* :class:`fwdpy11.genetic_values.SlocusGBR`

The above functions calculate the :math:`G` component of :math:`P`.  Custom trait value functions can be written in C++
in a similar manner as for custom fitness functions for standard "popgen" simulations.

Adding noise to trait values
----------------------------------------------------------


Mapping trait values to fitness
----------------------------------------------------------

In order to calulate fitness, we need a Python callable relating trait values to fitness.  In this section, I show the
implementation of two classes included in fwdpy11.

Stabilizing selection with constant parameters
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The first class models Gaussian stabilizing selection ("GSS") with respect to a constant optimum trait value and constant
intensity of selection. Under GSS models, :math:`w=e^{-\frac{(P-O)^2}{2VS}}`, where :math:`P` is the trait value,
:math:`O` is the optimum, and :math:`VS` determines the intensity of selection against deviations from the optimum
(*e.g* selection against variance in trait values). 

This class is rather simple:

* The constructor takes the optimum and :math:`VS` as parameters.
* The call function takes a diploid's trait value and uses the class data to return :math:`w=e^{-\frac{(P-O)^2}{2VS}}`.

Stabilizing selection with a moving optimum
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
    import fwdpy11.genetic_values
    import fwdpy11.wright_fisher
    import fwdpy11.gsl_random as gsl
    import fwdpy11.model_params
    import numpy as np

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

    p = {'nregions':[],
    'sregions':[fp11.GaussianS(0,1,1,0.25)],
    'recregions':[fp11.Region(0,1,1)],
    'rates':(0.0,2e-3,1e-3),
    'demography':np.array([N]*N,dtype=np.uint32),
    'gvalue':(fwdpy11.genetic_values.SlocusAdditive,(2.0,)),
    'gv2w':(fwdpy11.genetic_values.GSSmo,(gss_params,)),
    }

    params = fp11.model_params.ModelParams(**p)

    fwdpy11.wright_fisher.evolve(rng,pop,params)

