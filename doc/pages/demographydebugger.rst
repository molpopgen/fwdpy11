.. _demographydebugger:

Debugging Demographic models
=================================================

.. autoclass:: fwdpy11.DemographyDebugger
   :members:

.. ipython:: python

    import fwdpy11
    import fwdpy11.demographic_models.IM

    pop = fwdpy11.DiploidPopulation([100], 1.)

    events, t1, t2 = fwdpy11.demographic_models.IM.two_deme_IM(
        pop.N, 0.3, 0.2, (2, 3.33), [0.01, 0.1])

    d = fwdpy11.DemographyDebugger(pop, events)
    print(d.report)
