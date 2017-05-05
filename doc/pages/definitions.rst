.. _definitions:

Definitions
======================================================================

We need to start out by defining some terms that will be used throughout the documentation.  

Region vs locus
-----------------------------------------------------------

fwdpy11 allows mutation and recombination parameters to vary along a genomic segment.  Let's call the genomic segment a
**locus**.  Within a locus, positions with different mutation and recombination parameters are called **regions**.

Regions have the following properties:

* They are parameterized by objects in the :class:`fwdpy11.regions.Region` class hierarchy.
* They may overlap.  For example, a "coding region" may be modeled as two overlapping regions, one in which neutral
  mutations occur, and another in which selected mutations occur.  There could even be more than one "selected
  region"--one for deleterious mutations and another for beneficial mutations.

See :ref:`regions` for details on defining regions in fwdpy11.

A locus contains the following properties:

* They contain regions.
* Regions specify how recombination and mutation occur in a gamete.  See :ref:`gametes`.
* Loci cannot overlap, but recombination can occur between them.

A "single locus" simulation of a single deme means that an object of type :class:`fwdpy11.fwdpy11_types.SlocusPop` is
being evolved.  Within that locus, you may have as many regions as you want.  You may even think of them as different
"genes" (aka loci), but that is a different locus concept from what we are using here. See :ref:`slocuspop` for more
details on this type.

Evolving a multi-locus system in a single deme means evolving a :class:`fwdpy11.fwdpy11_types.MlocusPop` object.  See
:ref:`mlocuspop` for details.

.. _genetic_values:

Genetic values, fitness, etc.
-----------------------------------------------------------
