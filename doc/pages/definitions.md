(definitions)=

# Definitions

We need to start out by defining some terms that will be used throughout the documentation.

(simtypes)=

## Simulation types

When it comes to the effect of mutations on fitness, we have two types:

```{eval-rst}
.. glossary:: 

    direct

        The mutation's effect size(s) affects fitness directly.
        This is the :math:`s` of standard population genetic models.

    indirect

        The mutation's effect size(s) affect one or more traits.
        The trait vector/matrix of a diploid then determines fitness.
```

Thus, we have two types of simulation:

```{eval-rst}
.. glossary:: 

    "Standard population genetic simulation"

        A simulation where mutations directly affect fitness

    "Quantitative genetic simulation"

        A simulation where mutations affect trait (phenotype) values first,
        which in turn affects fitness.

```

## Genetic values, fitness, etc.

Ultimately, the simulations work by sampling parents proportional to their fitnesses
and generating offspring. Because fwdpy11 supports both "standard population genetic"
simulations as well as simulation of quantitative traits,
we need to be precise about the terms we're using.

We apply the following definitions here:

```{eval-rst}
.. glossary:: 

    genetic value
        A *genetic value*, or :math:`G`, is the result of applying a function
        to the gamete data at a locus.

    trait value
        A *trait value*, or *phenotype*, :math:`P`, reflects any adjustment
        made to :math:`G`.  For example, the addition of random noise such
        that :math:`P=G+E`, where :math:`E` is the noise.

    fitness
        A *fitness*, or :math:`w` results from applying a function mapping
        :math:`P` to fitness. :math:`w` is a non-negative float.
```

Consider the following two cases:

First,  standard population genetic models of the sort that [sfs_code][sfs_code] or [SLiM2][slim2]
are typically used for.  Mutations affecting fitness interact multiplicatively.
Here, mutations directly affect fitness.  In our terms, {math}`G = w`.

Next, we think about a quantitative trait under Gaussian stabilizing selection
with respect to an optimum.  Mutations interact additively (or multiplicatively) to
generate {math}`G` and the final trait value is {math}`P = G + N(0,\sigma)`,
reflecting environmental variation with mean zero and standard deviation
{math}`\sigma`.  For an optimum trait value {math}`O`, fitness is calculated as

```{math}

w = e^{-\frac{(O-P)^2}{2VS}},

```

where {math}`VS` reflects the intensity of selection against extreme values
of {math}`P`.

We can see from these two examples that some modeling scenarios allow us to go
straight from a diploid's data to fitness while others require multiple functions
to go from genotype to genetic value to trait value and then, finally, to
fitness.

More details on these topics can be found in:


    Give links to where details can be found.
    Discuss fitness having default value of 1 and traits 0.
```

## Stateful vs stateless genetic value calculations

A genetic value calculation that only requires a diploid, a gamete container,
and a mutation container as arguments is considered "stateless".
In contrast, if a calculation requires knowledge of the rest of the state of
the population, or somehow depends on an externally-defined object,
then it is "stateful".  For example, if fitness depends on the mean
genetic distance to all other individuals in the population,
then that is something that would need to be updated and
recorded each generation, making genetic value calculations "stateful".
Another example is the snowdrift model, which is shown in {ref}`pysnowdrift`.

[sfs_code]: http://sfscode.sourceforge.net/SFS_CODE/index/index.html

[slim2]: https://messerlab.org/slim/


