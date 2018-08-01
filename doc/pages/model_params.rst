.. _model_params:

Parameterizing a simulation
======================================================================

.. versionadded:: 0.1.1

.. versionchanged:: 0.2.0

    Reduced from many classes to a single class.

The parameters for a simulation are represented by the following classes:

* :class:`fwdpy11.model_params.ModelParams`

See the documentation for :mod:`fwdpy11.model_params` for more detail.

These classes contain getter/setter properties allowing the retrieval and assignment, respectively, of simulation
parameters.

The names of these properties are also the set of acceptable keyword arguments that can be passed to class constructors.

These classes also attempt to check that paremeters are properly defined.  Currently, such checking is incomplete and
will be improved in future versions.

.. note::
    We have made every attempt to allow instances of these objects
    to be pickled, so that they can be sent out to other processes.  
    However, they are only as pickle-able as the 
    types that they contain.  If you cannot pickle a ModelParams 
    object, send the raw parameters out to other processes and 
    then create the ModelParams objects.

In general, it will be easiest to construct instances of :class:`fwdpy11.model_params.ModelParams` with a :class:`dict`.

The following example sets up a standard population simulation with multiplicative fitness effects:

.. testcode::

    import fwdpy11
    import fwdpy11.model_params
    import fwdpy11.genetic_values
    import fwdpy11.wright_fisher
    import numpy as np

    popsize = 1e4
    theta = 10000
    rho = theta
    pdict = {'nregions': [fwdpy11.Region(0,1,1)],
                'sregions': [fwdpy11.ExpS(0,1,1,-1.0/(2.0*popsize),0.5)],
                'recregions': [fwdpy11.Region(0,1,1)],
                'rates': (theta/(4.*popsize),1e-4,rho/(4.*popsize)),
                'demography': np.array([popsize]*10, dtype = np.uint32),
                'gvalue': fwdpy11.genetic_values.SlocusMult(1.0)
            }
    params = fwdpy11.model_params.ModelParams(**pdict)
    params.validate()
    pop = fwdpy11.SlocusPop(int(popsize))
    rng = fwdpy11.GSLrng(42)
    fwdpy11.wright_fisher.evolve(rng, pop, params)

The key points are:

* You provide *lists* of neutral, selected, and recombination regions.  See :ref:`regions`.
* You provide a *tuple* of rates representing the **total** neutral mutation, selected mutation, and recombination
  rates, respectively.
* The genetic value function is specified by a *tuple*.  The first element of the tuple is a **type name** and the
  second element is the constructor arguments required to generate an instance of that type.  See the next example 
  for extending this idea to simulations of traits and random effects.

This example sets up a multi-locus simulation of a quantitative trait under Gaussian stabilizing selection with random
effects:

.. testcode::

    import fwdpy11
    import fwdpy11.model_params
    import fwdpy11.genetic_values
    import fwdpy11.genetic_value_noise
    import fwdpy11.wright_fisher
    import fwdpy11.multilocus
    import numpy as np
    import inspect
    popsize = 1e4
    theta = 10000
    rho = theta
    locus_boundaries = [(float(i),float(i)+1.0) for i in range(5)]

    gv2w = fwdpy11.genetic_values.GSS(VS=1.0,opt=0.0)
    noise = fwdpy11.genetic_value_noise.GaussianNoise(mean=0.0, sd=0.1)
    gvalue = fwdpy11.genetic_values.MlocusAdditive(1.0,gv2w,noise)
    pdict = {'nregions': [[fwdpy11.Region(i[0],i[1],1)] for i in locus_boundaries],
                'sregions': [[fwdpy11.GaussianS(i[0],i[1],1,0.1)] for i in locus_boundaries],
                'recregions': [[fwdpy11.Region(i[0],i[1],1)] for i in locus_boundaries],
                'rates': ([theta/(4.*popsize)]*len(locus_boundaries),
                        [1e-4]*len(locus_boundaries),
                        [rho/(4.*popsize)]*len(locus_boundaries)),
                'interlocus_rec': fwdpy11.multilocus.binomial_rec([0.5]*(len(locus_boundaries)-1)),
                'demography': np.array([popsize]*10, dtype = np.uint32),
                'prune_selected':False,
                'gvalue':gvalue
            }
    params = fwdpy11.model_params.ModelParams(**pdict)
    params.validate()
    pop = fwdpy11.MlocusPop(int(popsize), locus_boundaries)
    rng = fwdpy11.GSLrng(42)
    fwdpy11.wright_fisher.evolve(rng, pop, params)

The key differences from the single-locus example are:

* You pass in lists of lists of regions.  There is one list per locus, and each list specifies variation in rates for
  each locus separately.
* You pass in tuples of lists of rates.
* You have to specify how recombination occurs *between* loci.  See :py:mod:`fwdpy11.multilocus`.
* As with the genetic value, we pass in tuples of type names plust constructor argumetns to specify the genetic 
  value to fitness map ('gv2w') and the random effects on trait values ('noise').

Examples:

* :ref:`qtraits1`

