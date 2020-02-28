.. _data_types:

Data types related to simulated populations
======================================================================

Random number generator
-----------------------------------------------------------

.. autoclass:: fwdpy11.GSLrng

    .. autoattribute:: __init__


Base class for populations
----------------------------------------------------------

.. autoclass:: fwdpy11.PopulationBase

.. note:: 

    This base class is imperfect, and contains some attributes
    and functions that assume diploidy.  These issues will be
    fixed in future releases, largely by moving the diploid-specific
    stuff to the derived classes. To hopefully avoid confusion,
    the oddly-placed attributes are documented in 
    :class:`fwdpy11.DiploidPopulation`.

.. autoattribute:: fwdpy11.PopulationBase.haploid_genomes

    .. note::

        The haploid_genomes container contains extant and extinct
        genomes. For the latter, :attr:`fwdpy11.HaploidGenome.n`
        equals 0.

.. autoattribute:: fwdpy11.PopulationBase.mutations

.. autoattribute:: fwdpy11.PopulationBase.mcounts

.. autoattribute:: fwdpy11.PopulationBase.tables

.. autoattribute:: fwdpy11.PopulationBase.genetic_values

.. autoattribute:: fwdpy11.PopulationBase.ancient_sample_genetic_values

Populations of diploids
----------------------------------------------------------

.. autoclass:: fwdpy11.DiploidPopulation
   :members:
    
Diploid Genotypes
----------------------------------------

This class contains two integers pointing to the 
genomes making up the individual.  The genomes themselves
are represented by :class:`fwdpy11.HaploidGenome`.

.. autoclass:: fwdpy11.DiploidGenotype

    .. autoattribute:: first
    .. autoattribute:: second

The diploid genotypes are stored in the following container:

.. autoclass:: fwdpy11.DiploidVector

For all intents and purposes, this container behaves as a standard
Python list, but its actual representation is a contiguous array
that is handled on the C++ side.

This type supports the Python buffer protocol, meaning that it
can be cast to a numpy structured array without a copy, and the attribute
names are used as field names. 

Diploid metadata
------------------------------------------------------------

.. autoclass:: fwdpy11.DiploidMetadata
    :members:

Metadata are stored in the following container:

.. autoclass:: fwdpy11.DiploidMetadataVector

For all intents and purposes, this container behaves as a standard
Python list, but its actual representation is a contiguous array
that is handled on the C++ side.

This type supports the Python buffer protocol, meaning that it
can be cast to a numpy structured array without a copy, and the attribute
names are used as field names. 

HaploidGenomes
-----------------------------------------------------------

Class :class:`fwdpy11.HaploidGenome` describes a gamete:

.. autoclass:: fwdpy11.HaploidGenome

    .. autoattribute:: n
    .. autoattribute:: mutations
    .. autoattribute:: smutations

HaploidGenomes are stored in the following container:

.. autoclass:: fwdpy11.HaploidGenomeVector

For all intents and purposes, this container behaves as a standard
Python list, but its actual representation is a contiguous array
that is handled on the C++ side.
    
Mutations 
-----------------------------------------------------------

A mutation is described by :class:`fwdpy11.Mutation`:   

.. autoclass:: fwdpy11.Mutation

    The position of a mutation is a floating-point value:

    .. autoattribute:: pos

    For simulations with a single effect size/selection coefficient,
    the value is a float and held in the field:

    .. autoattribute:: s

    Likewise, the heterozygous effect/dominance of the variant is 
    in the attribute:

    .. autoattribute:: h

    For simulations with multivariate effects, the analogs of `s`
    and `h` are stored as numpy arrays:

    .. autoattribute:: esizes

    .. autoattribute:: heffects

    Other attributes of mutations include:

    .. autoattribute:: g

    .. autoattribute:: label

    .. note::

        The `label` attribute is assigned when a mutation is generated.

Mutations are stored in the following container:

.. autoclass:: fwdpy11.MutationVector

For all intents and purposes, this container behaves as a standard
Python list, but its actual representation is a contiguous array
that is handled on the C++ side.

Types related to tree sequence recording
==========================================================

.. autoclass:: fwdpy11.TableCollection
    :members:

.. autoclass:: fwdpy11.EdgeTable

.. autoclass:: fwdpy11.Edge

    .. autoattribute:: parent

    .. autoattribute:: child

    .. autoattribute:: left

    .. autoattribute:: right

.. autoclass:: fwdpy11.IndexedEdge

    .. autoattribute:: pos
    .. autoattribute:: time
    .. autoattribute:: parent
    .. autoattribute:: child

.. autoclass:: fwdpy11.NodeTable
.. autoclass:: fwdpy11.Node

    .. versionchanged:: 0.6.0

         The `population` field is renamed `deme`
    
    .. autoattribute:: deme
    .. autoattribute:: time

    
.. autoclass:: fwdpy11.MutationTable
.. autoclass:: fwdpy11.MutationRecord

    .. autoattribute:: key

    .. autoattribute:: node

    .. autoattribute:: site

    .. autoattribute:: derived_state

    .. autoattribute:: neutral

.. autoclass:: fwdpy11.SiteTable
.. autoclass:: fwdpy11.Site

    .. autoattribute:: position
    .. autoattribute:: ancestral_state

.. autoclass:: fwdpy11.TreeIterator
    :members:

.. autoclass:: fwdpy11.VariantIterator
    :members:

.. autoclass:: fwdpy11.DataMatrixIterator

    The class constructor is:
    
    .. autoattribute:: __init__

    The following properties are numpy arrays
    providing read-only access to an internal
    instance of :class:`fwdpy11.DataMatrix`:

    .. autoattribute:: neutral
    .. autoattribute:: selected
    .. autoattribute:: neutral_keys
    .. autoattribute:: selected_keys
    .. autoattribute:: neutral_positions
    .. autoattribute:: selected_positions

.. autoclass:: fwdpy11.NoAncientSamples

    .. autoattribute:: __init__

Types related to discrete demographic events
=======================================================

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

.. autoclass:: fwdpy11.SetDemeSize

   .. automethod:: __init__

   The following attributes give read-only access to instance data:

   .. autoattribute:: when
   .. autoattribute:: deme
   .. autoattribute:: new_size
   .. autoattribute:: resets_growth_rate

.. autoclass:: fwdpy11.SetSelfingRate

   .. automethod:: __init__

   The following attributes give read-only access to instance data:

   .. autoattribute:: when
   .. autoattribute:: deme
   .. autoattribute:: S

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

Miscellaneous types
=======================================================

.. autoclass:: fwdpy11.RecordNothing

    .. autoattribute:: __init__


