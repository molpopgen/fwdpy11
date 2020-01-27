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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The model here is one of soft selection [Levene1953]_, meaning that the number of 
breeding individuals ("adults") in each deme is fixed at a certain value.
A nice overview of this model and how it compares to others in [Felsenstein1976]_.
You may also find [Christiansen1974]_ and [Christiansen1975]_ useful.

Each generation, offspring ("juveniles") are generated in each deme.  Parents are drawn
from demes according to a migration matrix, if one is provided, else they are drawn from
the offspring deme.  Within a parental deme, a specific parent is chosen proportionally
to relative fitness within the deme.

The timings of events
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Below, we discuss various events that may happen.  These event types
include things like deme size changes, "mass migration" events, etc..
These events will occur at a certain time in a simulation. That time
refers to the birth time of a generation and the events are applied
*prior* to generating offspring, meaning that the events happen *to
the parents*.  For example, if half of deme zero moves and colonizes
a new deme (deme 1), then that means that half of the current alive individuals
(possible parents) have their ``deme`` field changed from zero to one
prior to generating any offspring.

.. _soft_sel_deme_setup:

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

Another method involves mass migration events at the beginning of a simulation.
See :ref:`massmigrations`.

The DiscreteDemography class
------------------------------------------------

The demographic events are stored in instances of :class:`fwdpy11.DiscreteDemography`.
These events, whose interface is described below, are passed in ``list`` objects
when created a :class:`fwdpy11.DiscreteDemography` instance.

These instances may be used to parameterize the ``demography`` field of a 
:class:``fwdpy11.ModelParams`` instance.  To illustrate this, here is a 
function that we'll use repeatedly below:


.. ipython:: python

    def setup_and_run_model(pop, ddemog, simlen, recorder=None, seed=654321):
        pdict = {'nregions': [],
                'sregions': [],
                'recregions': [],
                'rates': (0, 0, 0,),
                'gvalue': fwdpy11.Multiplicative(2.),
                'demography': ddemog,
                'simlen': simlen
               }
        params = fwdpy11.ModelParams(**pdict)
        rng = fwdpy11.GSLrng(654321)
        fwdpy11.evolvets(rng, pop, params, 100, recorder)


We will also define a simple class to record all deme sizes over time:


.. ipython:: python

    class SizeTracker(object):
        def __init__(self):
            self.data = []
        def __call__(self, pop, sampler):
            md = np.array(pop.diploid_metadata, copy=False)
            self.data.append((pop.generation, pop.N,
                             np.unique(md['deme'], return_counts=True)))


Event types
------------------------------------------------

The following sub-sections describe the various types of demographic
events allowed during a simulation.

.. _massmigrations:

Mass migrations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Mass migration events represent the "bulk" movement of individuals
in a single generation.  Such events allow you to model population
splits, merges, etc..

These events are represented by instances
of :class:`fwdpy11.MassMigration`.  Currently, you create instances
of this type using one of the following two functions:

* :func:`fwdpy11.copy_individuals`
* :func:`fwdpy11.move_individuals`

As the name implies, the first function creates an event that *copies*
individuals from a source deme to a destination.  The latter *moves*
them.

Both functions take five arguments, which may be used either named
or unnamed.  In order, they are:

* ``when``: the time (generation) when the event will occur
* ``source``: the ID of the source deme
* ``destination``: the ID of the destination deme
* ``fraction``: the fraction (proportion) of ``source`` moved/copied to ``dest``.
* ``resets_growth_rate``: If ``True``, the event resets the growth rate to :attr:`fwdpy11.NOGROWTH`
  in **both** ``source`` and ``dest``. If ``False``, growth rates remain unchanged.
  The default is ``False``.

These operations act on proportions of populations rather than on numbers
of individuals. Multiple events in a single generation are allowed, see
:ref:`multiple_mass_migrations`.

Setting the initial state of a simulation
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Let's look at an example where we use mass migration events to set up
"who is where" at the start of a simulation.  Since events happen in
the *parental* generation, we can use mass migrations to set up 
what demes individuals are in by applying events at generation 0.

The main difference between this method and that shown in
:ref:`soft_sel_deme_setup` is that these events move or copy *random*
individuals to new demes whereas using the  ``__init__`` approach 
builds the individuals in each deme sequentially.

For example, if we wish to start a simulation with 50 individuals in 
demes 0 and 50 in deme 1, we have two options:

1. Start with 50 individuals and *copy* them to deme 1 in generation 0
2. Start with 100 individuals and *move half of* them to deme 1 in generation 0

Here is the version implemented via a  copy:

.. ipython:: python

    pop = fwdpy11.DiploidPopulation(50, 1.)
    copy = [fwdpy11.copy_individuals(when=0, source=0, destination=1, fraction=1.0)]
    ddemog = fwdpy11.DiscreteDemography(mass_migrations=copy)
    setup_and_run_model(pop, ddemog, 1)
    md = np.array(pop.diploid_metadata, copy=False)
    np.unique(md['deme'], return_counts=True)


Here is what our object looks like:

.. ipython:: python

    print(copy[0])


Here is the version using a move:

.. ipython:: python

    pop = fwdpy11.DiploidPopulation(100, 1.)
    move = [fwdpy11.move_individuals(0, 0, 1, 0.5)]
    ddemog = fwdpy11.DiscreteDemography(mass_migrations=move)
    setup_and_run_model(pop, ddemog, 1)
    md = np.array(pop.diploid_metadata, copy=False)
    np.unique(md['deme'], return_counts=True)


For comparison, here is the object specifying the move:

.. ipython:: python

    print(move[0])

.. _multiple_mass_migrations:

Multiple mass migrations 
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

To specify multiple events, simply add more events to your list.
The events to not have to be sorted in any specific way.  Any sorting 
requirements get handled internally.

Multiple events involving the same source population in the same generation
need some explaining.   If the events are copies, things will tend to "just
work":

.. ipython:: python

    pop = fwdpy11.DiploidPopulation(50, 1.)
    copy = [fwdpy11.copy_individuals(0, 0, 1, 1.0),
            fwdpy11.copy_individuals(0, 0, 2, 1.0)]
    ddemog = fwdpy11.DiscreteDemography(mass_migrations=copy)
    setup_and_run_model(pop, ddemog, 1)
    md = np.array(pop.diploid_metadata, copy=False)
    np.unique(md['deme'], return_counts=True)
    

When the events are moves, it is not possible to move more than 100% 
of the individuals.  Attempting to do so will raise a ``ValueError``
exception:

.. ipython:: python

    pop = fwdpy11.DiploidPopulation(50, 1.)
    # Move all of deme 0 into demes 1 and 2,
    # which means we're trying to move 200% 
    # of deme 0...
    move = [fwdpy11.move_individuals(0, 0, 1, 1.0),
            fwdpy11.move_individuals(0, 0, 2, 1.0)]
    # ... which is not allowed
    try:
       ddemog = fwdpy11.DiscreteDemography(mass_migrations=move)
    except ValueError as e:
       print(e)

The rate of drift
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Moving versus copying individuals is an important modeling choice.
When you move individuals from one deme to another, the rate of drift
changes in the source deme (as its size is reduced).  This reduction
in size is also a sudden bottleneck.

Copying, on the other hand, does not change the rate of drift in the source
deme.  However, it does seem to imply some sudden increase in fecundity that
both came from nowhere and was short-lived.

The distinction between the two choices matters when trying to implement
models whose parameters have been inferred via approaches based on the
coalescent or on diffusion approximations. Such approaches usually assume
that population splits are copy events.

Instantaneous deme size changes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Instantaneous changes in deme size are managed by instances of 
:class:`fwdpy11.SetDemeSize`.

This class is relatively straightforward to use, so let's dive right in:

.. ipython:: python

    pop = fwdpy11.DiploidPopulation([20, 20], 1.)
    dd = fwdpy11.DiscreteDemography(set_deme_sizes=[fwdpy11.SetDemeSize(when=5,deme=1,new_size=100)])
    st = SizeTracker()
    setup_and_run_model(pop, dd, 10, st)
    for i in st.data:
        print(i)

You may also kill off demes by setting their size to zero:

.. ipython:: python

    pop = fwdpy11.DiploidPopulation([20, 20, 20], 1.)
    dd = fwdpy11.DiscreteDemography(set_deme_sizes=[fwdpy11.SetDemeSize(when=5,deme=1,new_size=0)])
    st = SizeTracker()
    setup_and_run_model(pop, dd, 6, st)
    for i in st.data:
        print(i)

Changing growth rates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Instances of :class:`fwdpy11.SetExponentialGrowth` manage the exponential growth rates per deme.
Growth rates less than one indicate population decline, greater than one means growth
and :attr:`fwdpy11.NOGROWTH` is equal to 1.0 to indicate no growth.

Changing the selfing rate
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


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

