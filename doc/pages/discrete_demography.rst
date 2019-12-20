.. _discretedemography:

.. ipython:: python
   :suppress:

   import fwdpy11

Discrete demography
======================================================================

.. autoclass:: fwdpy11.DiscreteDemography

   Representation of demographic events acting on 
   discrete demes (sub-populations).
   
   .. versionadded:: 0.5.3

   The constructor methods for this class are overloaded:

    .. automethod:: __init__

   The following attributes give read-only access to instance data:

    .. autoattribute:: mass_migrations   
    .. autoattribute:: set_growth_rates   
    .. autoattribute:: set_deme_sizes   
    .. autoattribute:: set_selfing_rates
    .. autoattribute:: migmatrix   
    .. autoattribute:: set_migration_rates   

All of the sections below should bring in the docstrings.

Changing growth rates
------------------------------------------------

.. data:: fwdpy11.NOGROWTH

   The author of this software has a bad habit of using the 
   value of 0.0 to represent "no exponential growth".  This is,
   of course, wrong, and so a constant is provided to help
   other people with similar problems:

.. ipython:: python

   print(fwdpy11.NOGROWTH)

.. autoclass:: fwdpy11.SetExponentialGrowth

   .. automethod:: __init__

   The following attributes give read-only access to instance data:

   .. autoattribute:: when
   .. autoattribute:: deme
   .. autoattribute:: G

* Pretty straightforward, but concrete examples will help a lot

Instantaneous deme size changes
------------------------------------------------

.. autoclass:: fwdpy11.SetDemeSize

   .. automethod:: __init__

   The following attributes give read-only access to instance data:

   .. autoattribute:: when
   .. autoattribute:: deme
   .. autoattribute:: new_size
   .. autoattribute:: resets_growth_rate

* How to determine if growth rate changes

Changing the selfing rate
------------------------------------------------

.. autoclass:: fwdpy11.SetSelfingRate

   .. automethod:: __init__

   The following attributes give read-only access to instance data:

   .. autoattribute:: when
   .. autoattribute:: deme
   .. autoattribute:: S

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

Mass migrations
------------------------------------------------

.. autoclass:: fwdpy11.MassMigration

   The following attributes give read-only access to instance data:

   .. autoattribute:: when
      
      When the event happens

   .. autoattribute:: source

      ID of the source deme

   .. autoattribute:: destination

      ID of the destination deme

   .. autoattribute:: fraction

      The fraction (proportion) of individuals migrating from `source`
      to `destination`.

   .. autoattribute:: move_individuals

      A boolean.  If `True`, then individuals are *moved* from `source`
      to destination.  If `False`, then they are copied.

   .. autoattribute:: resets_growth_rate

      A boolean.  If `True`, then the event resets the growth rate
      of *both* `source` and `destination` to :data:`fwdpy11.NOGROWTH`.


.. autofunction:: fwdpy11.move_individuals

.. autofunction:: fwdpy11.copy_individuals



* Move vs copy
* Proportions of demes
* Whether or not they change growth rates
* Merges and splits: examples
* This is how new demes are created
* Ghost population example


.. _migration:

Migration
------------------------------------------------

.. autoclass:: fwdpy11.MigrationMatrix

   .. automethod:: __init__

   The following attributes give read-only access to instance data:

   .. autoattribute:: shape

      The shape of the migration matrix

   .. autoattribute:: M

      Returns a copy of the rate matrix

   .. autoattribute scaled

      `True` if rates are multiplied by deme sizes during a simulation.

.. autoclass:: fwdpy11.SetMigrationRates

   .. automethod:: __init__

   The following attributes give read-only access to instance data:

   .. autoattribute:: when
   .. autoattribute:: deme
   .. autoattribute:: migrates

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
