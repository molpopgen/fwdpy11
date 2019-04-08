.. _genetic_values:


Genetic values
======================================================================

For some background, see :ref:`definitions`.

This section discusses the classes used to model genetic value and fitness calculations.  For the impatient, 
some examples probably suffice:

To simulate a "standard population genetic" simulation with additive effects on fitness:

.. testcode::

    import fwdpy11
    a = fwdpy11.Additive(2.0)
    print(a.is_fitness)
    print(a.scaling)

.. testoutput::

    True
    2.0

In the above, the `2.0` is the scaling of homozygous mutant fitness.  See the class documentation for
:class:`fwdpy11.Additive` for details.

To simulate fitness that is multiplicative across sites, use :class:`fwdpy11.Multiplicative`:

.. testcode::

    import fwdpy11
    m = fwdpy11.Multiplicative(2.0)
    print(m.is_fitness)
    print(m.scaling)


.. testoutput::

    True
    2.0

To simulate an additive-effects phenotype where phenotypes are subject to selection, we need to include
objects in the :class:`fwdpy11.GeneticValueToFitnessMap` class hierarchy.  For example, to model
Guassiant stablizing selection around a constant optimum trait value of zero with :math:`VS=1`, we use
:class:`fwdpy11.GSS`:

.. testcode::

    import fwdpy11
    a = fwdpy11.Additive(2.0, fwdpy11.GSS(0, 1))
    print(a.is_fitness)

.. testoutput::

    False

To model an optimum that changes over time, we can use :class:`fwdpy11.GSSmo`, which is used in :ref:`gss`.

.. autoclass:: fwdpy11.GeneticValue

.. autoclass:: fwdpy11.GeneticValueWithMapping

    Inherits from :class:`fwdpy11.GeneticValue`

.. autoclass:: fwdpy11.Additive
    :members:

    Inherits from :class:`fwdpy11.GeneticValueWithMapping`

.. autoclass:: fwdpy11.Multiplicative
    :members:

    Inherits from :class:`fwdpy11.GeneticValueWithMapping`

.. autoclass:: fwdpy11.GBR
    :members:

    Inherits from :class:`fwdpy11.GeneticValueWithMapping`

    .. testcode::

        import fwdpy11
        g = fwdpy11.GBR(fwdpy11.GSS(0, 1))
        g = fwdpy11.GBR(fwdpy11.GSS(0, 1), fwdpy11.GaussianNoise(0, 0.1))

.. autoclass:: fwdpy11.StrictAdditiveMultivariateEffects
    :members:

