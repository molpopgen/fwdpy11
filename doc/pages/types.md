---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

(data-types)=

# Data types related to simulated populations

## Random number generator

```{eval-rst}
.. autoclass:: fwdpy11.GSLrng

    .. autoattribute:: __init__

```

## Populations of diploids

```{eval-rst}
.. autoclass:: fwdpy11.DiploidPopulation
   :members:

   .. autoattribute:: fwdpy11.DiploidPopulation.generation
   .. autoattribute:: fwdpy11.DiploidPopulation.mutations
   .. autoattribute:: fwdpy11.DiploidPopulation.haploid_genomes
```

## Diploid Genotypes

This class contains two integers pointing to the
genomes making up the individual.  The genomes themselves
are represented by {class}`fwdpy11.HaploidGenome`.

```{eval-rst}
.. autoclass:: fwdpy11.DiploidGenotype

    .. autoattribute:: first
    .. autoattribute:: second
```

The diploid genotypes are stored in the following container:

```{eval-rst}
.. autoclass:: fwdpy11.DiploidVector
```

For all intents and purposes, this container behaves as a standard
Python list, but its actual representation is a contiguous array
that is handled on the C++ side.

This type supports the Python buffer protocol, meaning that it
can be cast to a numpy structured array without a copy, and the attribute
names are used as field names.

## Diploid metadata

```{eval-rst}
.. autoclass:: fwdpy11.DiploidMetadata
    :members:
```

Metadata are stored in the following container:

```{eval-rst}
.. autoclass:: fwdpy11.DiploidMetadataVector
```

For all intents and purposes, this container behaves as a standard
Python list, but its actual representation is a contiguous array
that is handled on the C++ side.

This type supports the Python buffer protocol, meaning that it
can be cast to a numpy structured array without a copy, and the attribute
names are used as field names.

## HaploidGenomes

Class {class}`fwdpy11.HaploidGenome` describes a gamete:

```{eval-rst}
.. autoclass:: fwdpy11.HaploidGenome

    .. autoattribute:: n
    .. autoattribute:: mutations
    .. autoattribute:: smutations
```

HaploidGenomes are stored in the following container:

```{eval-rst}
.. autoclass:: fwdpy11.HaploidGenomeVector
```

For all intents and purposes, this container behaves as a standard
Python list, but its actual representation is a contiguous array
that is handled on the C++ side.

## Mutations

A mutation is described by {class}`fwdpy11.Mutation`:

```{eval-rst}
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
```

Mutations are stored in the following container:

```{eval-rst}
.. autoclass:: fwdpy11.MutationVector
```

For all intents and purposes, this container behaves as a standard
Python list, but its actual representation is a contiguous array
that is handled on the C++ side.

# Types related to tree sequence recording

```{eval-rst}
.. autoclass:: fwdpy11.TableCollection
    :members:
```

```{eval-rst}
.. autoclass:: fwdpy11.EdgeTable
```

```{eval-rst}
.. autoclass:: fwdpy11.Edge

    .. autoattribute:: parent

    .. autoattribute:: child

    .. autoattribute:: left

    .. autoattribute:: right
```

```{eval-rst}
.. autoclass:: fwdpy11.NodeTable
```

```{eval-rst}
.. autoclass:: fwdpy11.Node

    .. versionchanged:: 0.6.0

         The `population` field is renamed `deme`

    .. autoattribute:: deme
    .. autoattribute:: time

```

```{eval-rst}
.. autoclass:: fwdpy11.MutationTable
```

```{eval-rst}
.. autoclass:: fwdpy11.MutationRecord

    .. autoattribute:: key

    .. autoattribute:: node

    .. autoattribute:: site

    .. autoattribute:: derived_state

    .. autoattribute:: neutral
```

```{eval-rst}
.. autoclass:: fwdpy11.SiteTable
```

```{eval-rst}
.. autoclass:: fwdpy11.Site

    .. autoattribute:: position
    .. autoattribute:: ancestral_state
```

```{eval-rst}
.. autoclass:: fwdpy11.TreeIterator
    :members:
```

```{eval-rst}
.. autoclass:: fwdpy11.VariantIterator
    :members:
```

```{eval-rst}
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
```

```{eval-rst}
.. autoclass:: fwdpy11.NoAncientSamples

    .. autoattribute:: __init__
```

# Types related to discrete demographic events

```{eval-rst}
.. autoclass:: fwdpy11.DiscreteDemography
    :members:
```

```{data} fwdpy11.NOGROWTH

The author of this software has a bad habit of using the
value of 0.0 to represent "no exponential growth".  This is,
of course, wrong, and so a constant is provided to help
other people with similar problems:

```

```{code-cell} python
import fwdpy11

print(fwdpy11.NOGROWTH)
```

```{eval-rst}
.. autoclass:: fwdpy11.SetExponentialGrowth
```

```{eval-rst}
.. autoclass:: fwdpy11.SetDemeSize
```

```{eval-rst}
.. autoclass:: fwdpy11.SetSelfingRate
```

```{eval-rst}
.. autoclass:: fwdpy11.MassMigration
```

```{eval-rst}
.. autofunction:: fwdpy11.move_individuals
```

```{eval-rst}
.. autofunction:: fwdpy11.copy_individuals
```

```{eval-rst}
.. autoclass:: fwdpy11.MigrationMatrix

   .. autoattribute:: shape

      The shape of the migration matrix

   .. autoattribute:: M

      Returns a copy of the rate matrix
```

```{eval-rst}
.. autoclass:: fwdpy11.SetMigrationRates
```

# Miscellaneous types

```{eval-rst}
.. autoclass:: fwdpy11.RecordNothing

    .. autoattribute:: __init__

.. autoclass:: fwdpy11.NewMutationData
   :members:

```


