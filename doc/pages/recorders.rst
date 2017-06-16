.. _recorders:

"Recorders" or "temporal samplers": taking time series data during a simulation
============================================================================================================================================

Background reading:

* :ref:`data_types`

Further reading once you're done here:

* :ref:`processingpops`
* :ref:`processingpopsNP`

fwdpy11 makes it easy to record data during a simulation.  Such recording is accomplished by callable objects that we
may call "recorders" or "temporal samplers".  fwdpy11 only implements one such object,
:class:`fwdpy11.temporal_samplers.RecordNothing`.  This object is used as a default in cases were no time series data
are needed.  You must define your own recorders, and they will all have the same common layout:

.. code-block:: python

    class MyRecorder(object):
        def __init__(self, **kwargs)
            """
            Class constructor takes parameters and 
            initializes class data.

            The use of the kwargs idiom here is purely for example.
            """
            pass
        def __call__(self, pop)
            """
            This will be called every generation and
            the population being simulated will be passed
            in
            """
            pass

In terms of details, that is about it.  You may interact with the population as described elsewhere in this manual,
updating the class data as appropriate.

I'll illustrate the concept with two complete examples.  The first example tracks population mean fitness over time.
We'll use numpy structured arrays to get mean fitness, as described in :ref:`processingpopsNP`.

.. code-block:: python

    import numpy as np

    class Wbar(object):
        def __init__(self):
            self.data=[]
        def __call__(self,pop):
            wbar = np.array(pop.diploids.trait_array())['w'].mean()
            self.data.append((pop.generation,wbar))

You may not want to record every generation.  You may want to record at a set of specific time points.  No problem:

.. code-block:: python

    import numpy as np

    class WbarAtTimes(object):
        def __init__(self,timepoints):
            self.timepoints=timepoints
            self.data=[]
        def __call__(self,pop):
            if len(self.timepoints) >0 and pop.generation == self.timepoints[0]:
                wbar = np.array(pop.diploids.trait_array())['w'].mean()
                self.data.append((pop.generation,wbar))
                self.timepoints.pop(0)

You may do *anything* with these objects involving valid Python data types and read-only access
to the populations.

Performance
-------------------------------------------------------------------------

Given that recorders will be called quite often, you want to make sure they are efficient.  See the following sections
of the manual:

* :ref:`processingpopsNP`
* :ref:`cython`
