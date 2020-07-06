.. _tskittransfer:

Data interchange with tskit
======================================================================

The tree sequence data structures may be converted to the analogous `tskit` objects using 
:func:`fwdpy11.DiploidPopulation.dump_tables_to_tskit`.

Let's take a look:

.. ipython:: python

   ts = pop.dump_tables_to_tskit()
   print(type(ts))


You may provide a `dict` that reflects the simulation parameters used.  This `dict`
will be part of the provenance information encoded in the :class:`tskit.TreeSequence`.
For example:

.. code-block:: python

    # Assuming mp is a fwdpy11.ModelParams:
    ts = pop.dump_tables_to_tskit(parameters={"model": str(mp), "seed": 12345})

Ultimately, it is up to you to decide what to include in `parameters`.
For example, if could be a script:

.. code-block:: python

    # Bonus points for somehow including the git commit hash corresponding
    # to the version of the script that you used!
    parameters = {"script": "/path/to/script", "type": "script", "seed": 1234}

In order to get the provenance information back out from a :class:`tskit.TreeSequence`:

.. code-block:: python

    import json

    provenance = json.decodes(ts.provenance(0).record)

If you recorded an instance of :class:`fwdpy11.ModelParams` as your `parameters`, you
can even reconstruct the original object (if you have the correct modules imported).
For example, if we assume that we encoded the model parameters as shown two listings ago:

.. code-block:: python

    import tskit
    import json
    import numpy as np
    import fwdpy11

    ts = tskit.load("sim.trees")
    provenance = json.decodes(ts.provenance(0).record)
    array = np.array # Annoyance!
    mp = eval(provenance["parameters"]["model"])

It is possible for model parameters to contain `numpy` arrays.  Unfortunately, their
string representations are not namespace-qualified, meaning that they say `array` rather
than `numpy.array` or `np.array`.
Thus, I made a type alias so that the `eval` would work.


Decoding the metadata
-----------------------------------------------

`tskit` and `fwdpy11` treat metadata quite differently.  The former is much more general, while the latter gives you direct access to the data objects on the C++ side.  The `tskit` approach is based on binary strings.  What `fwdpy11` does is encode strings that can be converted back to Python dictionaries.  For example, here is how one may process the individual metadata after dumping the tables to `tskit`:

.. ipython:: python

   import tskit

   individuals = ts.tables.individuals
   md = tskit.unpack_bytes(individuals.metadata, individuals.metadata_offset)
   first_individual = eval(md[0])
   print(type(first_individual))
   print(first_individual)

In order to distinguish "alive" from "dead" individuals (*e.g.*, those 
preserved as ancient samples), we need to make use of flags found in
:mod:`fwdpy11.tskit_tools`. For example, to identify all currently alive
individuals in the individual table:

.. ipython:: python

    import fwdpy11.tskit_tools

    alive_individuals = (
        ts.tables.individuals.flags & fwdpy11.tskit_tools.INDIVIDUAL_IS_ALIVE
    )

To find the nodes corresponding to those individuals:

.. ipython:: python

    alive_individual_nodes = (ts.tables.nodes.individual >= 0) & (
        (
            ts.tables.individuals.flags[ts.tables.nodes.individual]
            & fwdpy11.tskit_tools.INDIVIDUAL_IS_ALIVE
        )
        > 0
    )

    # The time of these nodes should all be zero because we simulated
    # non-overlapping generations
    alive_node_times = np.unique(
        ts.tables.nodes.time[np.where(alive_individual_nodes)], return_counts=True
    )
    alive_node_times

.. todo::

    Provide nice functions to return nodes from various times,
    etc..

The mutation metadata are trickier still, but the same general principles apply:

.. ipython:: python

   mutations = ts.tables.mutations
   sites = ts.tables.sites
   md = tskit.unpack_bytes(mutations.metadata, mutations.metadata_offset)
   print(eval(md[0]))
   # To get the position, we must use the sites table:
   print(sites[mutations[0].site].position)
   # Compare to what fwdpy11 has stored:
   m = pop.mutations[pop.tables.mutations[0].key]
   print(
       f"position: {m.pos}\ns: {m.s}\nh: {m.h}\n"
       f"label: {m.label}\n"
       f"g: {m.g}\n"
       f"neutral: {m.neutral}\n"
       f"esizes: {m.esizes}\n"
       f"heffects: {m.heffects}\n"
   )

Note that the `tskit` representation of the mutation record's the allele's
*age* in generations while :attr:`fwdpy11.Mutation.g` is the generation when
the mutation arose.  The reason for this discrepancy is because the two packages
record time in different directions!  The conversion to and from is trivial:

.. ipython:: python

   print(f"Origin to age = {pop.generation - m.g + 1}")

