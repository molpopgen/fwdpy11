.. _softselection:

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

This chapter describes how to generate detailed demographic models of individuals
evolving withing discrete sub-populations, or demes.

If you want to skip the details and see basic models provided by ``fwdpy11``, then
see :ref:`demographic_models`.

The model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The model here is soft selection [Levene1953]_, meaning that the number of 
breeding individuals ("adults") in each deme is fixed at a certain value.
A nice overview of this model and how it compares to others in [Felsenstein1976]_.
You may also find [Christiansen1974]_ and [Christiansen1975]_ useful.

Each generation, offspring ("juveniles") are generated in each deme.  Parents are drawn
from demes according to a migration matrix, if one is provided, else they are drawn from
the parents in the offspring deme.  Within a parental deme, a specific parent is
chosen proportionally to relative fitness within the deme.

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
    pop.deme_sizes()
    for m in pop.diploid_metadata:
       for n in m.nodes:
            assert m.deme == pop.tables.nodes[n].deme

Another method involves mass migration events at the beginning of a simulation.
See :ref:`massmigrations`.

The DiscreteDemography class
------------------------------------------------

The demographic events are stored in instances of :class:`fwdpy11.DiscreteDemography`.
These events, whose interface is described below, are passed in ``list`` objects
when creating a :class:`fwdpy11.DiscreteDemography` instance.

These instances may be used to parameterize the ``demography`` field of a 
:class:`fwdpy11.ModelParams` instance.  To illustrate this, here is a 
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
        rng = fwdpy11.GSLrng(seed)
        fwdpy11.evolvets(rng, pop, params, 100, recorder)


We will also define a simple class to record all deme sizes over time:


.. ipython:: python

    class SizeTracker(object):
        def __init__(self):
            self.data = []
        def __call__(self, pop, sampler):
            self.data.append((pop.generation, pop.N,
                             pop.deme_sizes()))


Compatibility with previous versions of fwdpy11
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Previous versions only supported size changes within a single deme.  These size changes were
parameterized via a ``numpy`` array specifying the size at each time point.  It is still possible
to specify the demography using that approach:

.. ipython:: python

       N = np.array([10]*10 + [5]*5 + [10]*10, dtype=np.uint32)
       pdict = {'nregions': [],
               'sregions': [],
               'recregions': [],
               'rates': (0, 0, 0,),
               'gvalue': fwdpy11.Multiplicative(2.),
               'demography': N
              }
       params = fwdpy11.ModelParams(**pdict)
       rng = fwdpy11.GSLrng(654321)
       pop = fwdpy11.DiploidPopulation(10, 1.0)
       fwdpy11.evolvets(rng, pop, params, 100)

Internally, the ``numpy`` array gets converted to instances of :class:`fwdpy11.SetDemeSize`, which is described
below (:ref:`set_deme_sizes`).  These instances are stored in a :class:`fwdpy11.DiscreteDemography` object:

.. ipython:: python

    print(params.demography)
    for i in params.demography.set_deme_sizes:
        print(i)

The simulation length is inferred from the ``numpy`` array, too:

.. ipython:: python

    params.simlen, len(N)

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

.. note::

   When a mass migration event *copies* individuals from deme, 
   the individuals copied are sampled *without replacement*.  Thus,
   if the fraction copied is 1.0, then every individual is copied.

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
    pop.deme_sizes()


Here is what our object looks like:

.. ipython:: python

    copy[0]


Here is the version using a move:

.. ipython:: python

    pop = fwdpy11.DiploidPopulation(100, 1.)
    move = [fwdpy11.move_individuals(0, 0, 1, 0.5)]
    ddemog = fwdpy11.DiscreteDemography(mass_migrations=move)
    setup_and_run_model(pop, ddemog, 1)
    pop.deme_sizes()


For comparison, here is the object specifying the move:

.. ipython:: python

    move[0]

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
    pop.deme_sizes()
    

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

.. _set_deme_sizes:

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

Let's look at an example:

.. ipython:: python

    pop = fwdpy11.DiploidPopulation([50], 1.)
    g = [fwdpy11.SetExponentialGrowth(when=0,deme=0,G=1.1)]
    dd = fwdpy11.DiscreteDemography(set_growth_rates=g)
    st = SizeTracker()
    setup_and_run_model(pop, dd, 6, st)
    for i in st.data:
        print(i)

The deme sizes each generation must be integer values.  The simulation uses C/C++ rules for
rounding double-precision values to integer values. The function ``numpy.rint`` uses the same
rules:

.. ipython:: python

   N0 = np.float(50.0)
   for i in range(6):
       Ni = N0*np.power(1.1,i+1)
       print(i+1, Ni, np.rint(Ni))

You may need to keep the rounding policy in mind when trying to predict final deme sizes when testing
or when trying to convert a model from continuous time into discrete time.

Changing the selfing rate
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Instances of :class:`fwdpy11.SetSelfingRate` affect the rate of selfing-versus-outcrossing in different
demes, or to change the rate within a deme over time. The default is that individuals don't self
unless they are picked twice as a parent by chance.

Using this type is straightforward.  Before we dive in, we will create a new recorder
type to track parents each generation:

.. ipython:: python

    class ParentTracker(object):
        def __init__(self):
            self.data = []
        def __call__(self, pop, sampler):
            for i in pop.diploid_metadata:
                 self.data.append((i.label, i.deme, i.parents))

Let's run a simulation for a couple of generations:
   
.. ipython:: python

    pop = fwdpy11.DiploidPopulation([5, 5], 1.)
    sr = [fwdpy11.SetSelfingRate(when=0, deme=1, S=1.0)] # Deme 1 always selfs
    dd = fwdpy11.DiscreteDemography(set_selfing_rates=sr)
    pt = ParentTracker()
    setup_and_run_model(pop, dd, 2, pt)

In our output, the deme label is the second value in each tuple, and any individual
in deme 1 has the same parent listed twice because they were the product of a selfing event:

.. ipython:: python

    for i in pt.data:
        print(i)

(In the above output, the parent IDs are the indexes of the parental individuals from their
generation.)

.. _migration:

Migration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For models with multiple demes, migration between then is managed by an
instance of :class:`fwdpy11.MigrationMatrix`.

For a migration matrix ``M``, the default interpretation of ``M[i, j]`` is the
fraction of deme ``i`` that will be replaced by migrations from deme ``j``. The 
entry ``M[i, i]`` represents the non-migrant fraction of deme ``i``'s ancestry.
The matrix is "row-major" meaning that rows refer to migration into source demes.
This definition of the migration matrix corresponds to that found in several
different sources ([Christiansen1974]_, [Christiansen1975]_).
This definition of migration is also what diffusion models assume (*e.g.* [Jouganous2017]_)
as well as coalescent simulations like *msprime* [Kelleher2016]_.

For example, consider the following matrix:

.. ipython:: python

   m = np.array([0.9, 0.1, 0.5, 0.5]).reshape(2,2)
   m

The first row corresponds to the ancestry of deme ``0``, such that 90% of parents will be
non-migrants and 10% will be migrants from deme ``1``:

.. ipython:: python

   m[0,]

To be concrete, if the size of deme ``0`` in the next generation is 1,000, then the expected
number of migrant and non-migrant parents of offspring in deme ``0`` is:

.. ipython:: python

   m[0,] * 1e3

The second row implies that half the ancestry of deme ``1`` is due to migrants and half
due to non-migrants:

.. ipython:: python

   m[1,]

The ``numpy`` array is sufficient to construct our demographic model:

.. ipython:: python

    d = fwdpy11.DiscreteDemography(migmatrix=m)
    d.migmatrix
    d.migmatrix.M

By default, there is no migration, which is represented by the value ``None``.  For example,
the following model has no migration events:

.. ipython:: python

    # Define demographic events w/o any migration stuff
    d = fwdpy11.DiscreteDemography(set_deme_sizes=[fwdpy11.SetDemeSize(0, 1, 500)])
    d.migmatrix is None

Likewise, if an identity matrix is provided an migration rates are never changed later,
then the input matrix is ignored:

.. ipython:: python

    d = fwdpy11.DiscreteDemography(migmatrix=np.identity(2))
    d.migmatrix is None

The only reason to use the identity matrix is to start a simulation with no migration
and then change the rates later via instances of :class:`fwdpy11.SetMigrationRates`.
To see this in action, we'll first generate a new type to track if parents of
offspring in deme 1 are migrants or not:

.. ipython:: python

    class MigrationTracker(object):
        def __init__(self, N0):
            self.N0 = N0
            self.data = []
        def __call__(self, pop, sampler):
            for i in pop.diploid_metadata:
                if i.deme == 1:
                    p = []
                    for j in i.parents:
                        if j < self.N0:
                            p.append((j, True))
                        else:
                            p.append((j, False))
                    self.data.append((pop.generation, i.label, i.deme, p))

.. ipython:: python

    # No migration at first
    mm = np.identity(2)
    # In generation 3, reset migration rates for deme 1 such
    # that parents are equally likey from both demes.
    cm = [fwdpy11.SetMigrationRates(3, 1, [0.5, 0.5])]
    dd = fwdpy11.DiscreteDemography(migmatrix=mm, set_migration_rates=cm)
    pop = fwdpy11.DiploidPopulation([10, 10], 1.0)
    mt = MigrationTracker(10)
    setup_and_run_model(pop, dd, 4, mt)

    for i in mt.data:
        nmig = 0
        if i[1] > 10:
            if i[3][0][1] is True:
                nmig+=1
            if i[3][1][1] is True:
                nmig+=1
        mstring = ""
        if nmig > 0:
            mstring="<- {} migrant parent".format(nmig)
        if nmig > 1:
            mstring += 's'
        print(i, mstring)

An alternative model of migration
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The description of migration rates above implies that migration events are 
independent of of source deme sizes.  To revisit our earlier example:

.. ipython:: python

   m = np.array([0.9, 0.1, 0.5, 0.5]).reshape(2,2)
   # The is the expected number of parents from demes 0 and 1
   # to offspring born in deme 0:
   m[0,] * 1000

``fwdpy11`` allows for a different migration scheme where the size of the source deme
matters.  For this model, ``M[i ,j]`` is the probability that an individual from
deme ``j`` is a parent in deme ``i``.  Internally, the migration matrix entries
``M[i, j]`` are multiplied by the size of the *source* demes, which means that
larger demes with nonzero migration rates to other demes have a larger chance
of being parents.

For example:

.. ipython:: python

   deme_sizes = np.array([1000, 2000])
   m
   md = m*deme_sizes
   # The following line divides each
   # row by its sum
   md/np.sum(md, axis=1)[:, None]

The first matrix is the same as in the preceding section--90% of the parents of deme
``0`` will be from deme ``0``.  In the second matrix, that fraction is reduced to
about 82% because deme ``1`` is twice as large as deme ``0``.

To enable this migration model, the following methods are equivalent:

.. ipython:: python
   
   # Method 1: pass a tuple with your numpy array and True
   # to indicate scaling M[i, j] by source deme sizes:
   d = fwdpy11.DiscreteDemography(migmatrix=(m, True))

   # Method 2: construct an instance of fwdpy11.MigrationMatrix,
   # passing True as the second argument to indicate the scaling
   # by source deme size.
   M = fwdpy11.MigrationMatrix(m, True)
   d = fwdpy11.DiscreteDemography(migmatrix=M)

.. note::

   This model of migration will typically give *different* results
   from diffusion models and coalescent simulations!


.. _migration_and_selfing:

Migration and selfing
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Within each deme, the selfing rate :math:`S` is the probability that an
individual selfs, and :math:`1-S` is the probability that an individual
outcrosses with another.

For a single deme, everything is very straightforward.  Likewise for many demes with no
migration.  The challenge arises when we have multiple demes, nonzero selfing rates in
one or more of them, and nonzero migration.

The challenge is due to the fact that  we consider the migration matrix elements
to be the probability of migration from deme ``j`` into deme ``i``.

If we focus on an offspring deme and pull a migrant parent from the migration matrix, one 
of two things may happen:

1. The migrant parent selfs, which occurs with probability :math:`S` for that migrant's deme.
2. The migrant parent outcrosses

In the second case, we have to go back to our migration matrix to choose another parent. Internally,
a second lookup table is used where each entry in row :math:`M_{i,-}` is multiplied by :math:`1 - S_j`,
where :math:`S_j` is the selfing probability in source deme :math:`j`.

Examples of models
-------------------------------------------------

Isolation with migration, or "IM"
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Consider two demes that split apart ``T`` time units ago and then grow to different
sizes in the present.  After the split, migration occurs between the two demes. The
demographic model has the following parameters:

* ``Nanc``, the ancestral population size.
* ``T``, the time of the split, which is in units of ``Nanc``.
* ``psplit``, the proportion of the ancestral population that splits off to found deme ``1``.
* ``N0``, the final size of deme ``0``, relative to ``Nanc``.
* ``N1``, the final size of deme ``1``, relative to ``Nanc``.
* ``m01``, the migration rate from deme ``0`` to deme ``1``.
* ``m10``, the migration rate from deme ``1`` to deme ``0``.

Here is the model in its entirety, with no mutation and no recombination.
First, we will set up the demographic events.  The population with evolve
for ``Nanc`` generations before the split.

.. ipython:: python

    Nanc = 100
    T = 0.2
    psplit = 0.33
    N0, N1 = 2, 3
    m01, m10 = 0.01, 0.0267

    # The split event
    split = [fwdpy11.move_individuals(when=Nanc, source=0,
                                      destination=1,
                                      fraction=psplit)] 
    # Get growth rates and set growth rate changes,
    # taking care to handle our rounding!
    gens_post_split = np.rint(Nanc*T).astype(int)
    N0split = np.rint(Nanc*(1.-psplit))
    N0final = np.rint(N0*Nanc)
    N1split = np.rint(Nanc*psplit)
    N1final = np.rint(N1*Nanc)
    G0 = fwdpy11.exponential_growth_rate(N0split, N0final,
                                         gens_post_split)
    G1 = fwdpy11.exponential_growth_rate(N1split, N1final,
                                         gens_post_split)
    growth = [fwdpy11.SetExponentialGrowth(Nanc, 0, G0),
              fwdpy11.SetExponentialGrowth(Nanc, 1, G1)]

    # Set up the migration matrix for two demes, but only 
    # deme zero exists.
    m = fwdpy11.migration_matrix_single_extant_deme(2, 0)
    # The rows of the matrix change at the split:
    cm = [fwdpy11.SetMigrationRates(Nanc, 0, [1.-m10, m10]),
          fwdpy11.SetMigrationRates(Nanc, 1, [m01, 1.-m01])]
    d = fwdpy11.DiscreteDemography(mass_migrations=split,
                                   set_growth_rates=growth,
                                   set_migration_rates=cm,
                                   migmatrix=m)

The above code made use of two helper functions:

* :func:`fwdpy11.exponential_growth_rate`
* :func:`fwdpy11.migration_matrix_single_extant_deme`

Finally, we can run it:

.. ipython:: python

    pop = fwdpy11.DiploidPopulation(Nanc, 1.0)
    setup_and_run_model(pop, d, Nanc+gens_post_split)

Now we check the final population sizes and make sure they are correct:

.. ipython:: python

    ds = pop.deme_sizes()
    assert ds[1][0] == N0final
    assert ds[1][1] == N1final

This model is common enough that you shouldn't have to implement it from 
scratch each time.  For this reason, we provide it in :func:`fwdpy11.demographic_models.IM.two_deme_IM`.

.. ipython:: python

    import fwdpy11.demographic_models.IM
    d2, tsplit, tafter_split = fwdpy11.demographic_models.IM.two_deme_IM(Nanc, T,
                                                                       psplit,
                                                                       (N0, N1),
                                                                       (m01, m10),
                                                                       burnin=1.0)
    pop2 = fwdpy11.DiploidPopulation(Nanc, 1.0)
    setup_and_run_model(pop2, d2, Nanc+gens_post_split)
    assert pop.generation == pop2.generation
    assert pop2.generation == tsplit + tafter_split
    ds2 = pop2.deme_sizes()
    assert np.array_equal(ds[0], ds2[0])
    assert np.array_equal(ds[1], ds2[1])

See :ref:`IMexample` for an example of using this function to compare to results
from diffusion models.

Run-time checking
-------------------------------------------------

The parameters of a demographic model are checked at run time at two different places:

* Upon object construction.  The various event objects try to make sure that the parameter inputs are valid.
* During a simulation. If invalid events occur during a simulation, the simulation raises a ``fwdpy11.DemographyError`` exception.

It is clearly preferable for a simulation to detect errors as early as possible.  While bad inputs can be
detected almost immediately, more subtle errors are only detected during simulation, which may take a while.
A more efficient approach to checking your models is described in :ref:`demographydebugger`.
