.. _data_types:

Data types related to simulated populations
======================================================================

This section describes the types that get updated when populations are evolved.  These are the same types that you 
access to analyze a simulated population. We will see specific examples of processing populations in :ref:`processingpops`.

Opaque Containers
-----------------------------------------------------------

Many data objects are stored in "opqaue" containers, which are thin wrappers around C++ containers.
These containers allow *direct* access the mutations allocated on the C++ side, meaning that there
is no copy from a C++ vector to a Python list, for example.

When an opaque container is "list-like", I will refer to it as an opaque list. Semantically, it is very similar to a
Python list, and you interact with it in the same way.  The main difference from a list is that negative indexing 
is not supported. However, they can be cast to a Python list if such operations are desired.

.. _popgenmuts:

Mutations 
-----------------------------------------------------------

A mutation is described by :class:`fwdpy11.Mutation`, which provides read-only access to the following
properties:

.. csv-table:: :class:`fwdpy11.Mutation` properties
    :header: "Property", "Definition"
    :widths: 5,5

    "pos", "Position"
    "s", "Selection coefficient/effect size"
    "h", "Dominance"
    "esizes", "A vector of effect sizes"
    "heffects", "A vector of heterozygous effects"
    "g", "Origination time"
    "neutral", "Boolean -- neutral or not?"
    "label", "16 bit integer"

Currently, all of the above fields are read-only except for "label".

Instances of :class:`fwdpy11.Mutation` are stored in an opaque list type called
:class:`fwdpy11.VecMutation`.  

.. note::

    A mutation container contains both *extinct* and *extant* mutations!  It is important to distinguish them.
    which is done using the mutation counts container (see next section).

    The reason why extinct mutations are present is because fwdpp (fwdpy11's C++ back-end) recycles the memory
    for extinct mutations in order to store new mutations.

.. _mcounts:

Mutation Counts
-----------------------------------------------------------

Mutation objects to not track their own frequencies.  Rather, they are stored in a
:class:`fwdpy11.VecUint32`, which is an opaque list of unsigned integers.

.. note::
    
    The length of a :class:`fwdpy11.VecUint32` in a simulation is the same
    length as the corresponding :class:`fwdpy11.VecMutation`.

.. note::

    Indices with values 0 zero correspond to the locations of extinct mutations in a mutation 
    container.

.. _mpos:

Mutation positions
-----------------------------------------------------------

.. versionadded:: 0.2.0

The positions of all extant variants, and their locations in the mutation container, are tracked via a list of tuples.
The first element in each tuple is the mutation position (a float) and the second it its location (an unsigned integer).

.. _gametes:

Gametes
-----------------------------------------------------------

Class :class:`fwdpy11.Gamete` describes a gamete.  It has the following read-only properties:

.. csv-table:: :class:`fwdpy11.Gamete` properties
    :header: "Property", "Definition"
    :widths: 5,5

    "n","Number of occurrences in population"
    "mutations","Container of keys to neutral mutations"
    "smutations","Container of keys to selected mutations"

The type of `mutations` and `smutations` is :class:`fwdpy11.VecUint32`, an opaque list of unsigned
integers.  These integers are the indexes of the mutations in the mutations container (and their counts in the mutation
counts container).

.. note::

    The `n` field does not imply that this precise gamete exists exactly `n` times in the population.  Rather, it refers
    to the number of times this specific instance exists.  The C++ back end does not require that unique gametes are
    represented once and only once.  If you want to know the frequency distribution at the level of gametes, you'd have
    to calculate that yourself by via an all-by-all comparison.

Gametes are stored in opaque lists of type :class:`fwdpy11.VecGamete`.

.. _diploids:

Diploids
-----------------------------------------------------------


.. versionchanged:: 0.2.0
    Meta data separated into a separate class from the gamete keys.

In a single-locus simulation, a diploid is represented by two classes, :class:`fwdpy11.DiploidGenotype`, and 
:class:`fwdpy11.DiploidMetadata`.  The former simply lists the two gametes making up a diploid:

.. csv-table:: :class:`fwdpy11.DiploidGenotype` properties
    :header: "Property", "Definition"
    :widths: 5,5

    "first", "Index of the first gamete."
    "second", "Index of the second gamete."

The later class contains useful additional data about the individual:


.. csv-table:: :class:`fwdpy11.DiploidMetadata` properties
    :header: "Property", "Definition"
    :widths: 5,5

    "w", "Fitness."
    "g", "Genetic value."
    "e", "Random component of trait value."
    "geography", "The x,y,z location of the diploid. (Not used in 0.2.0)"
    "label", "The index of this diploid in the population."
    "deme", "The deme index of this diploid. New in 0.2.0."
    "sex", "The sex index of this diploid. New in 0.2.0."
    "parents", "A list containing the label fields of each parent."
    "nodes", "The nodes in a tree sequence corresponding to this diploid. (Not used in 0.2.0)"

For a multi-locus simulation, the diploid genotype at each locus is stored in a :class:`fwdpy11.VecDiploid`, which is an opaque list of :class:`fwdpy11.DiploidGenotype` objects.  

In a single-locus simulation, diploids are stored in an opaque list of type
:class:`fwdpy11.VecDiploid`.  For multi-locus simulations, diploids are stored in
:class:`fwdpy11.VecVecDiploid`, which is also an opaque list.

.. _population:

The population base class
-----------------------------------------------------------

.. versionadded:: 0.2.0

.. versionchanged:: 0.3.0

    Added `genetic_values` and `ancient_samaple_genetic_values`.

All populations based around :class:`fwdpy11.Mutation` and :class:`fwdpy11.Gamete` inherit from a common base class,
:class:`fwdpy11.Population`.  This class in an Abstract Base Class, or ABC. You may not create instances of this class. \
Rather, you work with the derived classes :class:`fwdpy11.SlocusPop` and :class:`fwdpy11.MlocusPop`.

.. csv-table:: :class:`fwdpy11.Population` properties
    :header: "Property", "Definition"
    
    "N", "Current population size."
    "generation", "Current generation."
    "mutations", "A :class:`fwdpy11.VecMutation`. See :ref:`popgenmuts`."
    "mcounts", "See :ref:`mcounts`."
    "mut_lookup", "See :ref:`mpos`. The locations refer to the mutations container."
    "gametes", "A :class:`fwdpy11.VecGamete`.  See :ref:`gametes`."
    "fixations", "A :class:`fwdpy11.VecMutation` storing fixations. See :ref:`popgenmuts`."
    "fixation_times", "A :class:`fwdpy11.VecUint32` storing fixation times."
    "genetic_values", "A 2d numpy array of genetic values."
    "ancient_samaple_genetic_values", "A 2d numpy array of genetic values corresponding to ancient samples."

.. note::

    When simulating structured populations, `N` refers to the total number of individuals in the "meta-population",
    and specific deme data for individuals is obtained through :attr:`fwdpy11.DiploidMetadata.deme`.
    
.. _slocuspop:

Single-locus population objects
-----------------------------------------------------------

To simulate a single locus in a single deme, you use :class:`fwdpy11.SlocusPop`.  Instances of this
class are constructed with a population size:

.. testcode::

    import fwdpy11 as fp11
    pop = fp11.SlocusPop(10000)
    print(pop.N)
    print(pop.generation)

.. testoutput::

    10000
    0

These objects can be pickled. See :ref:`pickling_pops`.

This class contains the following read-only properties, in addition to those found in the base class
:class:`fwdpy11.Population`:

.. csv-table:: :class:`fwdpy11.SlocusPop` properties
    :header: "Property", "Definition"
    
    "diploids", "A :class:`fwdpy11.VecDiploid`.  See :ref:`diploids`."

.. _mlocuspop:

Multi-locus population objects
-----------------------------------------------------------

The type :class:`fwdpy11.MlocusPop` is analagous to :class:`fwdpy11.SlocusPop` in all but
one respect.  The `diploids` property type is :class:`fwdpy11.VecVecDiploid`.  See :ref:`diploids`
for details.  The class has all of the properties of :class:`fwdpy11.Population` plus the following:

.. csv-table:: :class:`fwdpy11.MlocusPop` properties
    :header: "Property", "Definition"

    "diploids", "A :class:`fwdpy11.VecVecDiploid`.  See :ref:`diploids`."
    "nloci", "The number of loci"
    "locus_boundaries", "The [begin,end) positions for each locus"

The need for `locus_boundaries` will be discussed elsewhere.

.. todo::

    Discuss locus boundaries somewhere.

Python data types stored in population objects
---------------------------------------------------------------------------------

.. versionadded:: 0.1.4

All population objects contain two generic Python objects.  These are called `popdata` and `popdata_user`.  The former
is read-only and exists to provide flexibility to the internal details of simulation functions.  The latter is a
read-write property that the user can modify.
