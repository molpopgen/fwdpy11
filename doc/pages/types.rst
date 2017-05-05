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

Mutation Containers
-----------------------------------------------------------

Instances of :class:`fwdpy11.fwdpp_types.Mutation` are stored in a container type called
:class:`fwdpy11.fwdpy11_type.MutationContainer`.  Semantically, this container behaves as if it were a Python list is
most respects.  The difference is that it *directly* access the mutations allocated on the C++ side, meaning that there
is no copy from a C++ vector to a Python list.  In pybind11_ terms, this is an "opaque" container, and fwdpy11 relies on
them for efficient access to the data types.

The main difference from a list is that negative indexing and slice indexing are not supported.

A mutation container can be cast to a list.


.. _pybind11: https://github.com/pybind/pybind11

