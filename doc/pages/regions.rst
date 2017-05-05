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

Mutation rates, recombination rates, and a weighting system
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A simulation will typically have a mutation rate, :math:`\\mu`, which represents the mean of a Poisson number of mutations per gamete per generation), and a recombination rate, :math:`r`, which again is the mean of Poisson number of crossover events (per diploid, per generation).  These parameters are the _total_ rates across an entire simulated region.  Variation in these parameters along the region are affected by a set of positions coupled with "weights", which the user specifies using S4 classes.

The base class: :class:`fwdpy11.regions.Region`

A :class:`fwdpy11.regions.Region` is a Python class with the following members:

* :math:`b`, which is the beginning/start of the region. The type is "float". 
* :math:`e`, which is the end/stop of the region. The type is "float".
* :math:`w`, which is a weighting factor associated with the region. The type is "float".

The members are used to inform the C++ code about the relative abundance of new mutations or recombination events will occur in what region.  Briefly, the number of events that occur in region :math:`i` are proportional to :math:`w_i/\sum_i w`, *i.e*, the weight assigned to region :math:`i` divided by the sum of weights assigned to all regions.  The weights for mutation events and for recombination events are considered separately.  Thus, in order to model a correlation between mutational processes and recombination, it is up to the user to generate regions whose weights are correlated.

fwdpy allows the :math:`w` slot to be interpreted in one of two ways:

* It is *not*  affected by the length of region.  Interally, the weight assigned is simply :math:`w`. 
* It is affected by the length of a region :math:`(e - b)`.

These two options are determined by arguments to class constructors, which we will see in examples below.  The latter is the default.

These two approaches allow for considerable modeling flexibility.  For example, the latter approach allows :math:`w` to be interpreted as a "per base-pair" rate.  Imagine that you wanted to simulate variation in recombination along discrete 100 kilobase chunks, and the rate of crossing-over *per base pair* increases in each chunk, and includes an initial chunk with no recombination:

1. start=1,stop= :math:`10^5`, :math:`r_{bp}=0`
2. start= :math:`10^5`,stop= :math:`2 \times 10^5`, :math:`r_{bp}=10^{-8}`
3. start= :math:`2 \times 10^5`,stop= :math:`3 \times 10^5`, :math:`r_{bp}=10^{-7}`  


This model boils down to the relative number of crossing overs per region occuring in the ratio :math:`0 : 10^{-8} : 10^{-7}`.  This is easily represented using fwdpy's classes:

.. testcode:: 

    import fwdpy11
    recRegions = [fwdpy11.Region(1,1e5,0),fwdpy11.Region(1e5,2e5,1e-8),fwdpy11.Region(2e5,3e5,1e-7)]
    for i in recRegions:
        print (i)


.. testoutput:: 

    beg = 1.000000000, end = 100000.000000000, weight = 0.000000000, label = 0
    beg = 100000.000000000, end = 200000.000000000, weight = 0.001000000, label = 0
    beg = 200000.000000000, end = 300000.000000000, weight = 0.010000000, label = 0


For this hypothetical example, the region lengths are all identical, and
thus an equivalent specification would be this:

.. testcode:: 

    recRegions = [fwdpy11.Region(1,1e5,0,False),fwdpy11.Region(1e5,2e5,1e-8,False),fwdpy11.Region(2e5,3e5,1e-7,False)]
    for i in recRegions:
        print (i)


.. testoutput::

    beg = 1.000000000, end = 100000.000000000, weight = 0.000000000, label = 0
    beg = 100000.000000000, end = 200000.000000000, weight = 0.000000010, label = 0
    beg = 200000.000000000, end = 300000.000000000, weight = 0.000000100, label = 0


Specific examples
-------------------

Mutations not affecting fitness ("neutral" mutations)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You specify regions where neutral mutations arise via the class :class:`fwdpy11.regions.Region`.  A region has a beginning, end, and a weight Thus, the following list would specify that 100% of neutral mutations occur on the continuous interval [0,1):

.. testcode::

    neutralRegions = [fwdpy11.Region(0,1,1)]

The beginning and end positions can be whatever you like:

.. testcode:: 

    #With a weight of 1, we're just rescaling the position here.
    neutralRegions = [fwdpy11.Region(0,100,1)]

To specify variation in the netural mutation process along a sequence,
combine multiple regions in your list:

.. testcode::

    #If coupled=False for the second region, the effect would be that region2's mutation rate per base pair is 10x less than region 1!!
    neutralRegions = [fwdpy11.Region(beg=0,end=1,weight=1),fwdpy11.Region(beg=2,end=12,weight=1,coupled=True)]

Internally, the total "mutational weight" of the first region will be a
function of its length, which is 1(1-0)=1. The second region's total
weight will be 1\*(12-2)=10, and it will have 10 times as many new mutations
arising as the first region.

.. testcode:: 

    #Let's see what happens if we set coupled=False:
    neutralRegions2 = [fwdpy11.Region(beg=0,end=1,weight=1),fwdpy11.Region(beg=2,end=12,weight=1,coupled=False)]
    print("The set with coupled=True:")
    for i in neutralRegions:
        print(i)
    print("The set with coupled=False:")
    for i in neutralRegions2:
        print(i)


.. testoutput::

    The set with coupled=True:
    beg = 0.000000000, end = 1.000000000, weight = 1.000000000, label = 0
    beg = 2.000000000, end = 12.000000000, weight = 10.000000000, label = 0
    The set with coupled=False:
    beg = 0.000000000, end = 1.000000000, weight = 1.000000000, label = 0
    beg = 2.000000000, end = 12.000000000, weight = 1.000000000, label = 0


See the difference in the above? (Look at the "weight" term in the
second line of each set.)

Mutations affecting fitness
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Type types of mutations affecting fitness that we consider will have two parameters associated with them:

* :math:`s`, the selection coefficient
* :math:`h`, the effect of the mutation in a heterozygote (a.k.a. the "dominance" of the mutation).

In a simulation, we may place a distribution on either :math:`s` itself or on the scaled selection parameter :math:`\alpha = 2Ns`.  These two methods are represented by the class :class:`fwdpy11.regions.Sregion`.  These classes contain/extend the :class:`fwdpy11.regions.Region` class described above, and thus inherit their members.  :class:`fwdpy11.regions.Sregion` adds :math:`h`, which is the dominance of a mutation, and then classes extending :class:`fwdpy11.regions.Sregion` add details about the distribution of fitness effects.  These classes are:

* :class:`fwdpy11.regions.ConstantS`
* :class:`fwdpy11.regions.UniformS`
* :class:`fwdpy11.regions.GammaS`
* :class:`fwdpy11.regions.GaussianS`
  
Crossover rate variation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Just like neutral mutations, intervals with different crossover rates are specified by different :class:`Region` objects.  Let's set up the following concrete example:

* A region where crossovers occur between positions [0,1)
* Positions [0,0.45) and [0.55,1) have uniform recombintion rates at the "background" rate.
* Positions [0.45,0.55) are a recombination hotspot with 100x the background intensity (per "base pair").

The above model can be represented as:

.. testcode::

    #recrate[2] is the hotspot:
    recrates = [fwdpy11.Region(0.,0.45,1.),fwdpy11.Region(0.55,1.,1.,),fwdpy11.Region(0.45,0.55,100.)]
    for i in recrates:
        print (i)


.. testoutput::

    beg = 0.000000000, end = 0.450000000, weight = 0.450000000, label = 0
    beg = 0.550000000, end = 1.000000000, weight = 0.450000000, label = 0
    beg = 0.450000000, end = 0.550000000, weight = 10.000000000, label = 0


Internally, this is what will happen to the above input:

* The total weight on the first region will be :math:`w = w \times (e-b) = 1\times(0.45-0) = 0.45`
* The weight on the second region will be :math:`1\times(1-0.55) = 0.45`
* The weight on the hotspot will be :math:`100\times(0.55-0.45) = 10`

This gives us what we want: the hotspot is 100x hotter "per base", and is 10% of the total region in length.  We therefore expect 10x as many crossovers in that region as in the flanking regions.

How to set up a model
---------------------------------

When setting up a model, it is important that you think in terms of conditional probabilities.  In other words, if the total rate to neutral variants is :math:`\mu_n`, then the weights passed along to a function have the interpretations "Given that a neutral mutation occurs, the probability that it occurs in a certain interval is :math:`x`", where :math:`x` is determined by the relative weight assigned to an interval.

The 'weights' that you assign are *relative* and need not sum to 1.  Each weight must be :math:`\geq 0`, though.

