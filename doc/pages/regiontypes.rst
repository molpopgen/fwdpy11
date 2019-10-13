Setting up genomic intervals
===================================================

.. autoclass:: fwdpy11.Region
    :members:

.. autoclass:: fwdpy11.Sregion
    
    This class is an ABC derived from :class:`fwdpy11.Region`

.. autoclass:: fwdpy11.ConstantS

    .. automethod:: __init__

.. autoclass:: fwdpy11.UniformS

    .. automethod:: __init__

.. autoclass:: fwdpy11.ExpS

    .. automethod:: __init__

.. autoclass:: fwdpy11.GammaS

    .. automethod:: __init__

.. autoclass:: fwdpy11.GaussianS

    .. automethod:: __init__

.. autoclass:: fwdpy11.MultivariateGaussianEffects

    .. automethod:: __init__

.. autoclass:: fwdpy11.GeneticMapUnit

.. autoclass:: fwdpy11.PoissonInterval

    This class inherits from :class:`fwdpy11.GeneticMapUnit`

    .. automethod:: __init__
    
.. autoclass:: fwdpy11.PoissonPoint

    This class inherits from :class:`fwdpy11.GeneticMapUnit`

    .. automethod:: __init__

.. autoclass:: fwdpy11.BinomialPoint

    This class inherits from :class:`fwdpy11.GeneticMapUnit`

    .. automethod:: __init__

.. autoclass:: fwdpy11.FixedCrossovers

    This class inherits from :class:`fwdpy11.GeneticMapUnit`

    .. automethod:: __init__

.. autoclass:: fwdpy11.BinomialInterval

    This class inherits from :class:`fwdpy11.GeneticMapUnit`

    .. automethod:: __init__
