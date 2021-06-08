(tskit_tools)=

# Module `fwdpy11.tskit_tools`

```{eval-rst}
.. automodule:: fwdpy11.tskit_tools
    :members:

    The following are boolean flags that are used to populate
    the ``flags`` colummn of a :class:`tskit.IndividualTable`.
    The flags represent individuals alive at the end of the simulation,
    those preserved as ancient samples, and those from the
    very first generation of the simulation, respectively.
    See :ref:`here <recapitation>` for why the last class
    of preserved samples are special.

    .. data:: fwdpy11.tskit_tools.INDIVIDUAL_IS_ALIVE
    .. data:: fwdpy11.tskit_tools.INDIVIDUAL_IS_PRESERVED
    .. data:: fwdpy11.tskit_tools.INDIVIDUAL_IS_FIRST_GENERATION

    .. autoclass:: fwdpy11.tskit_tools.WrappedTreeSequence

       This class has the following attributes:

       .. autoattribute:: fwdpy11.tskit_tools.WrappedTreeSequence.data

            Return the `data` field from top-level metadata

       .. autoattribute:: fwdpy11.tskit_tools.WrappedTreeSequence.demes_graph

            Return the `demes_graph` field from top-level metadata

       .. autoattribute:: fwdpy11.tskit_tools.WrappedTreeSequence.generation

            Return the `generation` field from top-level metadata

       .. autoattribute:: fwdpy11.tskit_tools.WrappedTreeSequence.model_params

            Return the `model_params` field from top-level metadata

       .. autoattribute:: fwdpy11.tskit_tools.WrappedTreeSequence.seed

            Return the `seed` field from top-level metadata

       .. autoattribute:: fwdpy11.tskit_tools.WrappedTreeSequence.ts

            Access the underlying :class:`tskit.TreeSequence`

       This class has the following methods:

       .. automethod:: fwdpy11.tskit_tools.WrappedTreeSequence.timepoints_with_individuals

       .. automethod:: fwdpy11.tskit_tools.WrappedTreeSequence.decode_individual_metadata

       .. automethod:: fwdpy11.tskit_tools.WrappedTreeSequence.decode_mutation_metadata

```

```{eval-rst}
.. autofunction:: fwdpy11.tskit_tools.decode_individual_metadata
```

```{eval-rst}
.. autoclass:: fwdpy11.tskit_tools.DiploidMetadata
```

```{eval-rst}
.. autofunction:: fwdpy11.tskit_tools.decode_mutation_metadata
```

