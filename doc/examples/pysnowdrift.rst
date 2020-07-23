.. _pysnowdrift:

Snowdrift model in python
===================================================================

A Python implementation of the model described in [Doebeli2004]_.
See :ref:`here <gvalues_python>` for how to implement genetic value
calculations in Python.

Implementing this model in C++ gives much better performance with little
difference in code complexity.  See :ref:`here <snowdrift_cpp_example>` for
the C++ version.

.. figure:: pysnowdrift.png

   Example output from the script.  The y axis values are simple
   indexes of sampling events. Sampling occurrs every 100 generations,
   so multiply by that value for the generation.

.. literalinclude:: ../../examples/python_genetic_values/pysnowdrift.py


