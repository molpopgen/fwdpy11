.. _tskittransfer:

Data interchange with tskit
======================================================================

The tree sequence data structures may be converted to the analogous `tskit` objects using the following function:

.. autofunction:: fwdpy11.DiploidPopulation.dump_tables_to_tskit

Let's take a look:

.. ipython:: python

   ts = pop.dump_tables_to_tskit()
   print(type(ts))

.. note::

   There is no function to dump a :class:`fwdpy11.TableCollection` directly.  The reason is that, if you
   have a new table collection that you got from simplifying within fwdpy11, you could have done the 
   exact same operation from within tskit. 

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

.. note::

   `tskit`'s IndividualTable does not distinguish "alive" from "dead" individuals.
   Thus, the metadata for :attr:`fwdpy11.DiploidPopulation.diploid_metadata` and
   :attr:`fwdpy11.DiploidPopulation.ancient_sample_metadata` are concatenated,
   with the former coming before the latter. For all nodes corresponding to individuals,
   their `flag` field is set to `tskit.NODE_IS_SAMPLE` in the nodes table.

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
   print(f"position: {m.pos}\ns: {m.s}\nh: {m.h}\n"
         f"label: {m.label}\n"
         f"g: {m.g}\n"
         f"neutral: {m.neutral}\n"
         f"esizes: {m.esizes}\n"
         f"heffects: {m.heffects}\n")

Note that the `tskit` representation of the mutation record's the allele's *age* in generations while :attr:`fwdpy11.Mutation.g` is the generation when the mutation arose.  The reason for this discrepancy is because the two packages record time in different directions!  The conversion to and from is trivial:

.. ipython:: python

   print(f"Origin to age = {pop.generation - m.g + 1}")
