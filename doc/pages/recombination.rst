Genetic maps
==============================================================

Recombination rates are allowed to vary along genomes in a discrete fashion.  fwdpy11 
provides two complementary methods for setting up such variation. 

By way of example, say we want the following genetic map:

* The total recombination rate per diploid is :math:`1e-3`, which is the mean of a Poisson process.
* Our genome is continuous on :math:`[0,10)`.
* The recombination rate is twice as high in :math:`[0, 5)` compared to :math:`[5, 10)`.

On approach uses an overall recombination rate and instances of :class:`fwdpy11.Region`:

.. testcode::

    import fwdpy11
    recrate = 1e-3
    recregions = [fwdpy11.Region(0, 5, 2),
                  fwdpy11.Region(5, 10, 1)]
    print(recregions[0])
    print(recregions[1])

.. testoutput::

    Region(beg=0, end=5, weight=10)
    Region(beg=5, end=10, weight=5)

The constructor arguments for the above are `start`, `stop`, and `weight`.

By default, a "region" object's total weight is recorded internally as 
:math:`weight\times (stop-start)`, which you can see in the output from the
print statements above.  Therefore it is helpful to view
:math:`weight` as the "rate per unit" and :math:`stop-start` as the 
number of units in the region.  This behavior can be changed, and the following
is equivalent to the above (because :math:`stop-start` is the same for each region):

.. testcode::

    recregions = [fwdpy11.Region(0, 5, 2, coupled=False),
                  fwdpy11.Region(5, 10, 1, coupled=False)]
    print(recregions[0])
    print(recregions[1])

.. testoutput::

    Region(beg=0, end=5, weight=2)
    Region(beg=5, end=10, weight=1)

Now, the `weight` arguments are treated as *absolute*, or exactly `2` and `1`.

In words, what we have is:

* The mean number of breakpoints per diploid is :math:`1e-3`
* For each breakpoint, its position is uniform on :math:`[0, 5)` with probability :math:`2/(2+1)`, or
  it is uniform on :math:`[5, 10)` with probability :math:`1/(2+1)`.

In essence, instances of :class:`fwdpy11.Region` parameterize a multinomial distribution that is used to 
choose the ranges within which breakpoints are uniformly-distributed.

The second approach uses instances of :class:`fwdpy11.PoissonInterval`:

.. testcode::
        
    import fwdpy11
    recregions = [fwdpy11.PoissonInterval(0, 5, 2e-3/3),
                  fwdpy11.PoissonInterval(5, 10, 1e-3/3)] 

The constructor arguments are `start`, `stop`, mean number of breakpoints.  Each breakpoint is uniform
on :math:`[start, stop)`.

Some comments are useful at this point:

* Different :math:`[start, stop)` intervals may overlap.  The interpretation of such a setup is your problem.
* The first method, based on :class:`fwdpy11.Region` is slightly faster, but less flexible.  More on the flexibility
  below.
* When using classes like :class:`fwdpy11.PoissonInterval`, the recombination rate that you use to construct a 
  :class:`fwdpy11.ModelParams` instance is ignored, as the rates are stored in the individual objects.

The class :class:`fwdpy11.PoissonInterval` inherits from the ABC :class:`fwdpy11.GeneticMapUnit`.  This class hierarchy
allows very flexible modeling of discrete variation in recombination rates.  The next example takes advantage of this
flexibility.

.. testcode::

    import fwdpy11
    NLOCI = 10
    LOCUS_LENGTH = 10
    RECRATE_PER_LOCUS = 1e-3
    LOCUS_BOUNDARIES = [(i,i+LOCUS_LENGTH) for i in range(0, NLOCI*LOCUS_LENGTH, LOCUS_LENGTH)]
    recregions = [fwdpy11.PoissonInterval(*i, RECRATE_PER_LOCUS) for i in LOCUS_BOUNDARIES]
    for i in LOCUS_BOUNDARIES[:-1]:
        recregions.append(fwdpy11.BinomialPoint(i[1], 0.5))


This example sets up 10 genomic regions, each 10 "units" long.  Within each region, the mean number of breakpoints (per 
diploid, per generation) is :math:`1e-3`.  Between each region, a single recombination occurs with probability of
one-half, meaning that each region is assorting independently (50cM between each region).

See :class:`fwdpy11.GeneticMapUnit` to learn about other classes in this hierarchy.

As an aside, this example is not creating objects in order by their positions.  Such ordering is not required.

.. note::

    In general, you probably do not want to specify "gaps" between regions when setting up genetic maps!
    If you do that, you are implying that there is a non-recombining region between two regions,
    and mutating the tables (see :func:`fwdpy11.infinite_sites`) will place neutral variants in the gaps!
