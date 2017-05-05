Data types
======================================================================

Mutations
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

Opaque Containers
-----------------------------------------------------------

Many data objects are stored in "opqaue" containers, which are thin wrappers around C++ containers.
These containers allow *direct* access the mutations allocated on the C++ side, meaning that there
is no copy from a C++ vector to a Python list, for example.

The opaque containers discussed here differ from a list in that negative indexing and slice indexing are not supported.
However, they can be cast to a Python list if such operations are desired.

When an opaque container is "list-like", I will refer to it as an opaque list.


Mutation Containers
-----------------------------------------------------------

Instances of :class:`fwdpy11.fwdpp_types.Mutation` are stored in an opaque list type called
:class:`fwdpy11.fwdpy11_types.MutationContainer`.  

.. note::

    A mutation container contains both *extinct* and *extant* mutations!  It is important to distinguish them.
    which is done using the mutation counts container (see next section).

    The reason why extinct mutations are present is because fwdpp (fwdpy11's C++ back-end) recycles the memory
    for extinct mutations in order to store new mutations.

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

Gamete containers
-----------------------------------------------------------

Gametes are stored in opaque lists of type :class:`fwdpy11.fwdpy11_types.GameteContainer`.

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
