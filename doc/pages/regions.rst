.. _regions:

Regions: modeling discrete variation in mutation and recombination
======================================================================

Please see :ref:`definitions` before reading this section.

The simulation routines allow the details of the mutation and recombination models to vary along a "sequence" or "region".  A user is able to specify the details of such variation by passing *lists* to package functions.  For example, you are able to:

* Vary the neutral mutation rate along a sequence.
* Vary the distribution of selection coefficients (and the dominance associated with selected mutations) along a sequence.
* Vary the recombination rate along a sequence.

The implementation of such variation along a region is *discrete*.  A region is specified by a beginning, and end, and a weight, plus any additional data required to specify selection coefficients, dominance, etc.

Background
--------------------------------------------------
The models are parameterized through Python's "new-style" class system.

Regions and weights
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A "region" refers to a continuous genomic interval associated with a weight.  The base class for a region is
:class:`fwdpy11.Region`, which contains the following members:

* :math:`b`, which is the beginning/start of the region. The type is "float". 
* :math:`e`, which is the end/stop of the region. The type is "float".
* :math:`w`, which is a weighting factor associated with the region. The type is "float".

The members are used to inform the C++ code about the relative abundance of new mutations or recombination events will occur in what region.  Briefly, the number of events that occur in region :math:`i` are proportional to :math:`w_i/\sum_i w`, *i.e*, the weight assigned to region :math:`i` divided by the sum of weights assigned to all regions.  The weights for mutation events and for recombination events are considered separately.  Thus, in order to model a correlation between mutational processes and recombination, it is up to the user to generate regions whose weights are correlated.

fwdpy11 allows the :math:`w` slot to be interpreted in one of two ways:

* It is *not*  affected by the length of region.  Interally, the weight assigned is simply :math:`w`. 
* It is affected by the length of a region :math:`(e - b)`.

These two options are determined by arguments to class constructors, which we will see in examples below.  The latter is the default.

.. note:: 
 
    The 'weights' that you assign are *relative* and need not sum to 1.  Each weight must be :math:`\geq 0`, though.

Variation in mutation rates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Mutations not affecting fitness ("neutral" mutations)
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. note:: 
    
   This section only applies to simulations *without* tree sequence recording.
   When simulating with tree sequence recording, neutral mutations are added *after*
   the simulation is complete.

You specify regions where neutral mutations arise via the class :class:`fwdpy11.Region`.  A region has a beginning, end, and a weight Thus, the following list would specify that 100% of neutral mutations occur on the continuous interval [0,1):

.. ipython:: python
    :suppress:
    
    import fwdpy11

.. ipython:: python

    neutralRegions = [fwdpy11.Region(0,1,1)]

The beginning and end positions can be whatever you like:

.. ipython:: python 

    #With a weight of 1, we're just rescaling the position here.
    neutralRegions = [fwdpy11.Region(0,100,1)]

To specify variation in the netural mutation process along a sequence,
combine multiple regions in your list:

.. ipython:: python

    #If coupled=False for the second region, the effect would be that region2's mutation rate per base pair is 10x less than region 1!!
    neutralRegions = [fwdpy11.Region(beg=0,end=1,weight=1),fwdpy11.Region(beg=2,end=12,weight=1,coupled=True)]

Internally, the total "mutational weight" of the first region will be a
function of its length, which is 1(1-0)=1. The second region's total
weight will be 1\*(12-2)=10, and it will have 10 times as many new mutations
arising as the first region.

.. ipython:: python

    #Let's see what happens if we set coupled=False:
    neutralRegions2 = [fwdpy11.Region(beg=0,end=1,weight=1),fwdpy11.Region(beg=2,end=12,weight=1,coupled=False)]
    print("The set with coupled=True:")
    for i in neutralRegions:
        print(i)
    print("The set with coupled=False:")
    for i in neutralRegions2:
        print(i)

See the difference in the above? (Look at the "weight" term in the
second line of each set.)

Mutations affecting fitness
++++++++++++++++++++++++++++++++

Type types of mutations affecting fitness that we consider will have two parameters associated with them:

* :math:`s`, the selection coefficient
* :math:`h`, the effect of the mutation in a heterozygote (a.k.a. the "dominance" of the mutation).

In a simulation, we may place a distribution on either :math:`s` itself or on the scaled selection parameter :math:`\alpha = 2Ns`.  These two methods are represented by the class :class:`fwdpy11.Sregion`.  These classes contain/extend the :class:`fwdpy11.Region` class described above, and thus inherit their members.  :class:`fwdpy11.Sregion` adds :math:`h`, which is the dominance of a mutation, and then classes extending :class:`fwdpy11.Sregion` add details about the distribution of fitness effects.  These classes are:

* :class:`fwdpy11.ConstantS`
* :class:`fwdpy11.UniformS`
* :class:`fwdpy11.GammaS`
* :class:`fwdpy11.ExpS`
* :class:`fwdpy11.GaussianS`

.. versionchanged:: 0.13.a2
    Added ability to have these DFE objects represent distributions of scaled selection parameter via the "scaling"
    attribute.
  
Variation in recombination rates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. versionchanged:: 0.3.0

    Update to discuss more general approach to genetic maps

There are two approaches to modeling variation in recombination rates.  

First, the simulation may specify a recombination rate, `r`, and then use a list of :class:`fwdpy11.Region`
to model variation along the genome. This is essentially the same approach described above for neutral mutations.

.. _generalized_maps: 

A more general approach to genetic maps
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. versionadded:: 0.3.0

An alternative approach to modeling variation in recombination rates involves classes derived from
:class:`fwdpy11.GeneticMapUnit` (which is an ABC).  These classes allow you to "compose" a genetic map that is a mixture of continuous
intervals and point processes.

The relevant classes are:

* :class:`fwdpy11.PoissonInterval`, which specifies that the number of breakpoints are Poisson-distributed and positions
  uniform on the continuous interval :math:`[beg, end)`.
* :class:`fwdpy11.FixedCrossovers` generates a fixed number of breakpoints whose positions are 
  uniform on the continuous interval :math:`[beg, end)`.
* :class:`fwdpy11.BinomialPoint` represents recombination events occurring at a specific position with a specific
  probability.
* :class:`fwdpy11.PoissonPoint` also represents recombination events occurring at a fixed position.  The number of
  breakpoints is Poisson-distributed, and a breakpoint is inserted if the total number is odd.

For example, to model two continuous regions separated by 25 centiMorgans:

.. ipython:: python

    recRegions = [fwdpy11.PoissonInterval(0, 1, 1e-3),
                  fwdpy11.BinomialPoint(1, 0.25),
                  fwdpy11.PoissonInterval(1, 2, 1e-3)]

To model a genomic segment having exactly one crossover on the interval :math:`[0,1)`:

.. ipython:: python

    recRegions = [fwdpy11.FixedCrossovers(0, 1, 1)]

The number of recombination breakpoints in the intervals :math:`[0,1)` and :math:`[1,2)` will both be Poisson-distributed with means
of :math:`10^{-3}`.  A recombination event *between* the two regions will happen in 25% of meioses.

.. note::

    This scheme differs from the above (based on Regions) in that there is no overall recombination rate that needs to be specified.
    Rather, the rates are used to construct the individual objects.


