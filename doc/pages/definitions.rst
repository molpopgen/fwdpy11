.. _definitions:

Definitions
======================================================================

We need to start out by defining some terms that will be used throughout the documentation.  

Region vs locus
-----------------------------------------------------------

fwdpy11 allows mutation and recombination parameters to vary along a genomic segment.  We apply the following definitions:

.. glossary:: 

    locus
        A genomic segment

    region
        A subset of a locus.  Recombination and mutation parameters may vary across regions. Regions may overlap.

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

Ultimately, the simulations work by sampling parents proportional to their fitnesses and generating offspring. Because
fwdpy11 supports both "standard population genetic" simulations as well as simulation of quantitative traits, we need to
be precise about the terms we're using.

We apply the following definitions here:

.. glossary::

    genetic value
        A *genetic value*, or :math:`G`, is the result of applying a function to the gamete data at a locus. 

    trait value
        A *trait value*, or *phenotype*, :math:`P`, reflects any adjustment made to :math:`G`.  For example, the addition of
        random noise such that :math:`P=G+E`, where :math:`E` is the noise. 

    fitness
        A *fitness*, or :math:`w` results from applying a function mapping :math:`P` to fitness. :math:`w` is a non-negative
        float.

Consider the following two cases:

First,  standard population genetic model of the sort that sfs_code_ or SLiM2_  are typically used for.  Mutations
affecting fitness interact multiplicatively.  Here, mutations directly affect fitness.  In our terms, :math:`G = w`.

Next, we thing about a quantitative trait under Gaussian stabilizing selection with respect to an optimum.  Mutations
interact multiplicatively to generate :math:`G` and the final trait value is :math:`P = G + N(0,\sigma)`, reflecting
environmental variation with mean zero and standard deviation :math:`\sigma`.  For an optimum trait value :math:`O`,
fitness is calculated as 

.. math::

    w = e^{-\frac{(O-P)^2}{2VS}},

where :math:`VS` reflects the intensity of selection against extreme values of :math:`P`. (See :ref:`heritability` for more 
on :math:`VS`)

We can see from these two examples that some modeling scenarios allow us to go straight from a diploid's data to fitness
while others require multiple functions to go from genotype to genetic value to trait value and then, finally, to
fitness.

More details on these topics can be found in:

* :ref:`model_params`
* :ref:`qtraits1`

Stateful vs stateless genetic value calculations
-----------------------------------------------------------

A genetic value calculation that only requires a diploid, a gamete container, and a mutation container as argumetnts is
considered "stateless".  In contrast, if a calculation requires knowledge of the rest of the state of the population, or
somehow depends on an externally-defined object, then it is "stateful".  For example, if fitness depends on the mean
genetic distance to all other individuals in the population, then that is something that would need to be updated and
recorded each generation, making genetic value calculations "stateful".  Another example is the snowdrift model, which
is shown in :ref:`stateful_fitness`.

.. _sfs_code: http://sfscode.sourceforge.net/SFS_CODE/index/index.html
.. _SLiM2: https://messerlab.org/slim/
