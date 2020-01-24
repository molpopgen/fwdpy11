.. softselection:

.. ipython:: python
   :suppress:

   import fwdpy11
   import numpy as np

Soft selection with discrete demes
======================================================================

This page describes implementing models of demographic events affecting
multiple demes.  This functionality was first released in version 0.6.0
and makes use of low-level types added in 0.5.3.

.. note::

   As of 0.6.0, these features only apply to simulations using tree sequence
   recording.

Overview
------------------------------------------------

The model
+++++++++++++++++++++++++++++++

The model here is one of soft selection [Levene1953]_, meaning that the number of 
breeding individuals ("adults") in each deme is fixed at a certain value.  A nice overview of this model and how it compares to others in [Felsenstein1976]_.  You
may also find [Christiansen1974]_ and [Christiansen1975]_ useful.

Each generation, offspring ("juveniles") are generated in each deme.  Parents are drawn
from demes according to a migration matrix, if one is provided, else they are drawn from the offspring deme.  Within a parental deme, a specific parent is chosen proportionally to relative fitness within the deme.

The timings of events
++++++++++++++++++++++++++++++

Below, we discuss various events that may happen.  These event types
include things like deme size changes, "mass migration" events, etc..
These events will occur at a certain time in a simulation. That time
refers to the birth time of a generation and the events are applied
*prior* to generating offspring, meaning that the events happen *to
the parents*.  For example, if half of deme zero moves and colonizes
a new deme (deme 1), then that means that half of the current alive individuals
(possible parents) have their ``deme`` field changed from zero to one
prior to generating any offspring.

Setting the initial demes in a simulation
------------------------------------------------

At the start of a simulation, you may assign diploids to demes 
when constructing an instance of :class:`fwdpy11.DiploidPopulation`.
For example, to initialize a population with 25 individuals in demes ``0`` and ``1``:

.. ipython:: python

    pop = fwdpy11.DiploidPopulation([25, 25], 1.0)
    md = np.array(pop.diploid_metadata, copy=False)
    np.unique(md['deme'], return_counts=True)
    print(np.unique([i.deme for i in pop.tables.nodes], return_counts=True))
    for m in pop.diploid_metadata:
       for n in m.nodes:
            assert m.deme == pop.tables.nodes[n].deme


Event types
------------------------------------------------

The following sub-sections describe the various types of demographic
events allowed during a simulation.

Mass migrations
++++++++++++++++++++++++++++++++++++++++++++++++


* Move vs copy
* Proportions of demes
* Whether or not they change growth rates
* Merges and splits: examples
* This is how new demes are created
* Ghost population example

Instantaneous deme size changes
++++++++++++++++++++++++++++++++++++++++++++++++

* How to determine if growth rate changes

Changing growth rates
++++++++++++++++++++++++++++++++++++++++++++++++


* Pretty straightforward, but concrete examples will help a lot


Changing the selfing rate
++++++++++++++++++++++++++++++++++++++++++++++++


* Straightforward.  Main thing to point out is that selfing competes with migration.  Still
  working out the mechanics here.

Current thoughts for how to deal w/selfing follow.

Within each deme, the selfing rate :math:`S` is the probability that an individual selfs,
and :math:`1-S` is the probability that an individual outcrosses with another.

For a single deme, everything is very straightforward.  Likewise for many demes with no
migration.  The challenge arises when we have multiple demes, nonzero selfing rates in
one or more of them, and nonzero migration.

The challenge is due to the fact that, right now, we consider the migration matrix elements
to be the probability of migration from deme `r` into deme `c`, multiplied by the current
size of deme `r`. Here, `r` and `c` mean `row` and `column`.

If we focus on an offspring deme and pull a migrant parent from the migration matrix, one 
of two things may happen:

1. The migrant parent selfs, which occurs with probability :math:`S` for that migrant's deme.
2. The migrant parent outcrosses

In the first case, we are done and we use the migrant parent "twice" to generate the offspring.

In the second case, we have to go back to our migration matrix, and now we have our problem.
With respect to our offspring deme, the relevant column in the migration matrix is the weighted rates
of migration from all demes into the offspring deme.  What we really need is the *probability of
an outcrossing event* being a parent in our offspring deme.  Thus, it seems we need a **second**
lookup table where the "raw" migration weights are all weighted by the current :math:`1-S` for
each source deme.

.. _migration:

Migration
++++++++++++++++++++++++++++++++++++++++++++++++


* The MigrationMatrix class
* The m-by-m migration matrix represents the probability that an offspring in column c has a parent from
  row r, and the matrix is consulted for each parent (barring selfing, see above).  Thus, rows are
  source demes, and columns are destination demes.  I think we can say that this is the same forward
  migration matrix as in Christiansen and others, 1970s, but will have to check.
* scaled vs "un-scaled"
* In general, the various ways the values can be interpreted: probabilities, rates, scaled or not
* Show the various ways to construct things/initiate via the DiscreteDemography class. None 
  means no migration.

Debugging Demographic models
-------------------------------------------------

TBD -- probably a later PR

