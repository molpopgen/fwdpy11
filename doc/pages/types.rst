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

Mutations with a single effect size
-----------------------------------------------------------

A mutation is described by :class:`fwdpy11.fwdpp_types.Mutation`, which provides read-only access to the following
properties:

.. csv-table:: :class:`fwdpy11.fwdpp_types.Mutation` properties
    :header: "Property", "Definition"
    :widths: 5,5

    "pos", "Position"
    "s", "Selection coefficient/effect size"
    "h", "Dominance"
    "g", "Origination time"
    "neutral", "Boolean -- neutral or not?"
    "xtra", "16 bit integer"

Instances of :class:`fwdpy11.fwdpp_types.Mutation` are stored in an opaque list type called
:class:`fwdpy11.fwdpy11_types.MutationContainer`.  

.. note::

    A mutation container contains both *extinct* and *extant* mutations!  It is important to distinguish them.
    which is done using the mutation counts container (see next section).

    The reason why extinct mutations are present is because fwdpp (fwdpy11's C++ back-end) recycles the memory
    for extinct mutations in order to store new mutations.

.. _mcounts:

Mutation Counts
-----------------------------------------------------------

Mutation objects to not track their own frequencies.  Rather, they are stored in a
:class:`fwdpy11.fwdpy11_types.VectorUint32`, which is an opaque list of unsigned integers.

.. note::
    
    The length of a :class:`fwdpy11.fwdpy11_types.VectorUint32` in a simulation is the same
    length as the corresponding :class:`fwdpy11.fwdpy11_types.MutationContainer`.

.. note::

    Indices with values 0 zero correspond to the locations of extinct mutations in a mutation 
    container.

.. _gametes:

Gametes
-----------------------------------------------------------

Class :class:`fwdpy11.fwdpy11_types.Gamete` describes a gamete.  It has the following read-only properties:

.. csv-table:: :class:`fwdpy11.fwdpp_types.Gamete` properties
    :header: "Property", "Definition"
    :widths: 5,5

    "n","Number of occurrences in population"
    "mutations","Container of keys to neutral mutations"
    "smutations","Container of keys to selected mutations"

The type of `mutations` and `smutations` is :class:`fwdpy11.fwdpy11_types.VectorUint32`, an opaque list of unsigned
integers.  These integers are the indexes of the mutations in the mutations container (and their counts in the mutation
counts container).

.. note::

    The `n` field does not imply that this precise gamete exists exactly `n` times in the population.  Rather, it refers
    to the number of times this specific instance exists.  The C++ back end does not require that unique gametes are
    represented once and only once.  If you want to know the frequency distribution at the level of gametes, you'd have
    to calculate that yourself by via an all-by-all comparison.

Gametes are stored in opaque lists of type :class:`fwdpy11.fwdpy11_types.GameteContainer`.

.. _diploids:

Diploids
-----------------------------------------------------------

In a single-locus simulation, a diploid is represented by :class:`fwdpy11.fwdpy11_types.SingleLocusDiploid`, which
contains the following read-only properties:

.. csv-table:: :class:`fwdpy11.fwdpy11_types.SingleLocusDiploid` properties
    :header: "Property", "Definition"
    :widths: 5,5

    "first", "Index of the first gamete."
    "second", "Index of the second gamete."
    "w", "Fitness."
    "g", "Genetic value."
    "e", "Random component of trait value."
    "label", "The index of this diploid in the population."

For a multi-locus simulation, the diploid genotype at each locus is stored in a :class:`fwdpy11.fwdpy11_types.DiploidContainer`, which is an opaque list of :class:`fwdpy11.fwdpy11_types.SingleLocusDiploid` objects.  **The w/g/e/label fields are only populated for the first locus.**

.. note::

    Future changes to fwdpp will likely make the storage of data in a multi-locus diploid more efficient and sensible.

In a single-locus simulation, diploids are stored in an opaque list of type
:class:`fwdpy11.fwdpy11_types.DiploidContainer`.  For multi-locus simulations, diploids are stored in
:class:`fwdpy11.fwdpy11_types.VecDiploidContainer`, which is also an opaque list.

.. _slocuspop:

Single-locus, single-deme population objects
-----------------------------------------------------------

To simulate a single locus in a single deme, you use :class:`fwdpy11.fwdpy11_types.SlocusPop`.  Instances of this
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

This class contains the following read-only properties:

.. csv-table:: :class:`fwdpy11.fwdpp_types.SlocusPop` properties
    :header: "Property", "Definition"
    
    "N", "Current population size."
    "generation", "Current generation."
    "mutations", "A :class:`fwdpy11.fwdpy11_types.MutationContainer`. See :ref:`popgenmuts`."
    "mcounts", "See :ref:`mcounts`."
    "gametes", "A :class:`fwdpy11.fwdpy11_types.GameteContainer`.  See :ref:`gametes`."
    "diploids", "A :class:`fwdpy11.fwdpy11_types.DiploidContainer`.  See :ref:`diploids`."
    "fixations", "A :class:`fwdpy11.fwdpy11_types.MutationContainer` storing fixations. See :ref:`popgenmuts`."
    "fixation_times", "A :class:`fwdpy11.fwdpp_types.VectorUint32` storing fixation times."

.. _mlocuspop:

Multi-locus, single-deme population objects
-----------------------------------------------------------

The type :class:`fwdpy11.fwdpy11_types.MlocusPop` is analagous to :class:`fwdpy11.fwdpy11_types.SlocusPop` in all but
one respect.  The `diploids` property type is :class:`fwdpy11.fwdpy11_types.VecDiploidContainer`.  See :ref:`diploids`
for details.  The class has all of the properties of :class:`fwdpy11.fwdpy11_types.SlocusPop` plus the following:

.. csv-table:: :class:`fwdpy11.fwdpp_types.MlocusPop` properties
    :header: "Property", "Definition"

    "nloci", "The number of loci"
    "locus_boundaries", "The [begin,end) positions for each locus"

The need for `locus_boundaries` will be discussed elsewhere.

.. todo::

    Discuss locus boundaries somewhere.
