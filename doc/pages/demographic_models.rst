.. _demographic_models:

The demography debugger
======================================================================

.. autoclass:: fwdpy11.DemographyDebugger
   :members:

Pre-computed demographic models
======================================================================

.. automodule:: fwdpy11.demographic_models
   :members:

.. autoclass:: fwdpy11.demographic_models.DemographicModelDetails

.. autoclass:: fwdpy11.demographic_models.DemographicModelCitation

``fwdpy11.demographic_models.IM``
----------------------------------------------------------------------

.. automodule:: fwdpy11.demographic_models.IM

.. autofunction:: fwdpy11.demographic_models.IM.two_deme_IM

.. autoclass:: fwdpy11.demographic_models.IM.TwoDemeIMParameters

.. autoclass:: fwdpy11.demographic_models.IM.TwoDemeIMMetaData

Human demographic models
----------------------------------------------------------------------

This module provides the following models:

* [Tennessen2012]_ via :func:`fwdpy11.demographic_models.human.tennessen`.
* The three deme model from Table 2 of [Jouganous2017]_ via :func:`fwdpy11.demographic_models.human.jouganous_three_deme`.

.. automodule:: fwdpy11.demographic_models.human
   :members:
