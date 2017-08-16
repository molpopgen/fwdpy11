.. _model_params:

Parameterizing a simulation.
======================================================================

.. versionadded:: 0.1.1

The parameters for a simulation are represented by the following classes:

* :class:`fwdpy11.model_params.ModelParams`
* :class:`fwdpy11.model_params.SlocusParams` 
* :class:`fwdpy11.model_params.SlocusParamsQ` 
* :class:`fwdpy11.model_params.MlocusParams` 
* :class:`fwdpy11.model_params.MlocusParamsQ` 

See the documentation for :mod:`fwdpy11.model_params` for more detail.

These classes contain getter/setter properties allowing the retrieval and assignment, respectively, of simulation
parameters.

The names of these properties are also the set of acceptable keyword arguments that can be passed to class constructors.

These classes also attempt to check that paremeters are properly defined.  Currently, such checking is incomplete and
will be improved in future versions.

.. note::
    We have made every attempt to allow instances of these objects
    to be pickled, so that they can be sent out to other processes.  
    However, they are only as pickle-able as the 
    types that they contain.  If you cannot pickle a ModelParams 
    object, send the raw parameters out to other processes and 
    then create the ModelParams objects.


Examples:

* :ref:`qtraits1`

