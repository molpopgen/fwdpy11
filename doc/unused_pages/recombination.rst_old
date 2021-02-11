.. _geneticmaps:

Genetic maps
==============================================================

Recombination rates are allowed to vary along genomes in a discrete fashion.  `fwdpy11` 
provides two complementary methods for setting up such variation. 

Method 1: regions and weights
-------------------------------------------------------------

The method described in this section works in combination with a total overall recombination
rate.  This rate is the mean of a Poisson distribution and the intervals where recombination
breakpoints happen are chosen based on their relative weights.  The regions are instances
of :class:`fwdpy11.Region`.  A region represents a continuous, half-open interval within which
crossover positions are uniformly distributed.

By way of example, say we want the following genetic map:

* The total recombination rate per diploid is :math:`1e-3`, which is the mean of a Poisson process.
* Our genome is continuous on :math:`[0,10)`.
* The recombination rate is twice as high in one part of the "genome" than in the other.

To initialize a region object, the following parameters may be used:

* ``beg``, the start of the region
* ``end``, the end of the region
* ``weight``, the weight assigned to the region
* ``coupled`` is a ``bool`` and determines how the weights are handled internally (see below).
* ``label`` is an integer that defaults to ``0`` and is not relevant to recombination.

The first three parameters are required.  A valid region has :math:`beg \geq 0`,
:math:`end > beg` and :math:`weight >= 0` and defines a half-open interval :math:`[beg, end)`.

.. note::

    A ``weight`` of ``0`` is the same as simply not defining a region! There is 
    no requirement that all genetic map elements cover the entire genome. We allow
    zero-weight regions for those who think that it is cleaner/more explicit to write
    them down.

By default, ``coupled=True``, which means that the *total* weight assigned to a region
will be :math:`weight\times (end-beg)`.  It is helpful to view :math:`weight` as
the "rate per unit" and :math:`end-beg` as the number of units in the region. (For example,
"unit" could refer to base pairs, but it need not.)

There are two ways to set this model up.  The first is arguably the most intuitive, which is to make
one region twice as long as the other:

.. ipython:: python

    import fwdpy11

    recrate = 1e-3
    recregions = [
        fwdpy11.Region(beg=0., end=10./3., weight=1),
        fwdpy11.Region(beg=10./3., end=10, weight=1),
    ]
    for r in recregions:
        print(r)

.. note::

     The `recrate` value above would be set for a simulation
     via the ``rates`` ``kwarg`` to :class:`fwdpy11.ModelParams`

In the output, you see that ``coupled=True``, which means that the simulation's
back-end will assign twice as many crossovers to the second region as to the first.

In words, the ``recregions`` list and ``recrate`` value mean the following:

* The number of crossovers per diploid is Poisson distributed with mean
  ``recrate``, or ``0.001``.
* Each crossover breakpoint has a ``1/3`` chance of being uniformly 
  distributed in :math:`[0, 10/3)` and a ``2/3`` chance of being
  uniformly distributed in :math:`[10/3, 10)`.

A more abstract approach relies on setting ``coupled=False``, which means
that the "raw" weights that you input are the exact values used internally:

.. ipython:: python

    recregions = [
        fwdpy11.Region(beg=0, end=5, weight=1, coupled=False),
        fwdpy11.Region(beg=5, end=10, weight=2, coupled=False),
    ]
    for r in recregions:
        print(r)

Now, the `weight` arguments are treated as *absolute*, or exactly `1` and `2`, respectively.

In words, what we have is:

* The number of breakpoints per diploid is Poisson distributed with mean :math:`1e-3`
* For each breakpoint, its position is uniform on :math:`[0, 5)` with probability :math:`2/(2+1)`, or
  it is uniform on :math:`[5, 10)` with probability :math:`1/(2+1)`.

In essence, instances of :class:`fwdpy11.Region` parameterize a multinomial distribution that is used to 
choose the ranges within which breakpoints are uniformly-distributed.  A limitation of this approach
is that we cannot model discrete jumps in genetic maps, such as those between chromosomes.

Method 2: using "genetic map" classes
---------------------------------------------------------------------------

.. versionadded:: 0.3.0

An alternate approach uses instances of classes derived from the `ABC`
:class:`fwdpy11.GeneticMapUnit`. Here `Unit` refers to an *element* of
a genetic map rather than the actual units (`cM`, etc.).  Instances of
these classes contain their own rates and we can mix and match regions
where recombination breakpoints are Poisson and binomially distributed.

Let's revisit the example from the previous section.  This time, we will
use :class:`fwdpy11.PoissonInterval`:

.. ipython:: python

    recregions = [
        fwdpy11.PoissonInterval(beg=0, end=5, mean=2e-3 / 3),
        fwdpy11.PoissonInterval(beg=5, end=10, mean=1e-3 / 3),
    ]

The number of breakpoints in each :math:`[beg, end)` interval is Poisson distributed
with the given mean. The position of each breakpoint is uniform on :math:`[beg, end)`.

These classes also allow us to specify breakpoints at a specific position with a specific probability.
The next example sets up 4 genomic regions, each 10 "units" long.  Within each region, the mean number of breakpoints (per 
diploid, per generation) is :math:`1e-3`.  Between each region, a single recombination occurs with probability of
one-half, meaning that each region is assorting independently (50 `cM` between each region).

.. ipython:: python

    NLOCI = 4
    LOCUS_LENGTH = 10
    RECRATE_PER_LOCUS = 1e-3
    LOCUS_BOUNDARIES = [
        (i, i + LOCUS_LENGTH) for i in range(0, NLOCI * LOCUS_LENGTH, LOCUS_LENGTH)
    ]
    recregions = [fwdpy11.PoissonInterval(*i, RECRATE_PER_LOCUS) for i in LOCUS_BOUNDARIES]
    for i in LOCUS_BOUNDARIES[:-1]:
        recregions.append(fwdpy11.BinomialPoint(i[1], 0.5))
    for i in recregions:
        print(i)

As an aside, this example is not creating objects in order by their positions.  Such ordering is not required.

The following classes are available:

* :class:`fwdpy11.PoissonInterval`
* :class:`fwdpy11.PoissonPoint`
* :class:`fwdpy11.BinomialInterval`
* :class:`fwdpy11.BinomialPoint`
* :class:`fwdpy11.FixedCrossovers`

General comments
-------------------------------------------------------------

* Different :math:`[beg, end)` intervals may overlap.  The interpretation of such a setup is your problem.
* The first method, based on :class:`fwdpy11.Region` is slightly faster, but less flexible.  More on the flexibility
  below.
* When using classes like :class:`fwdpy11.PoissonInterval`, the recombination rate that you use to construct a 
  :class:`fwdpy11.ModelParams` instance is ignored, as the rates are stored in the individual objects.
* You do not need to specify regions with zero recombination. Their existence is implied given the total
  length of the genome being simulated (:attr:`fwdpy11.TableCollection.genome_length`).

.. note::

    Adding neutral mutations to the tables with :func:`fwdpy11.infinite_sites` will place
    neutral variants in the non-recombining regions.

