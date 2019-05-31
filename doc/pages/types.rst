.. _data_types:

Data types related to simulated populations
======================================================================

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

Populations of diploids
----------------------------------------------------------

.. autoclass:: fwdpy11.DiploidPopulation
    
    This class inherits from :class:`fwdpy11.PopulationBase`.

.. autoattribute:: fwdpy11.DiploidPopulation.diploids

    This field is a C++ container containing instances
    of :class:`fwdpy11.DiploidGenotype`.

Individuals are associated with metadata, which is represented by :class:`fwdpy11.DiploidMetadata`.
The metadata are stored separately for the current generation and for ancient samples:

.. autoattribute:: fwdpy11.DiploidPopulation.diploid_metadata

.. autoattribute:: fwdpy11.DiploidPopulation.ancient_sample_metadata



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

    .. autoattribute:: edges
    .. autoattribute:: nodes
    .. autoattribute:: mutations
    .. autoattribute:: genome_length


.. autoclass:: fwdpy11.EdgeTable

.. autoclass:: fwdpy11.Edge

    .. autoattribute:: parent

    .. autoattribute:: child

    .. autoattribute:: left

    .. autoattribute:: right

.. autoclass:: fwdpy11.NodeTable
.. autoclass:: fwdpy11.Node
    
    .. autoattribute:: population

    .. autoattribute:: time

    
.. autoclass:: fwdpy11.MutationTable
.. autoclass:: fwdpy11.MutationRecord

    .. autoattribute:: key

    .. autoattribute:: node

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
