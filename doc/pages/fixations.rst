.. _handling_fixations:

Handling fixations during a simulation
===========================================================================================

Background reading:

* :ref:`data_types`
* :ref:`simtypes`
* :ref:`genetic_values`
* :ref:`genetic_values_types`

This section discusses how fixations (mutations at a relative frequency of 1) are treated during simulations.  This
section is simple, and it also matters quite a bit!  The gist of this section is that we currently have a tradeoff between modeling the biology
and the efficiency of a simulation.

The biology is as follows:

* For standard population genetic simulations, relative fitness matters, and relative fitness is constant *up to a
  multiplicative constant*
* For simulations of traits, the absolute distance of a trait value from the optimum determines fitness.

The software engineering issue is:

* The model parameters objects have a *prune_selected* field, which is a boolean.  When `True`, mutations affecting
  fitness are recorded in a populations `fixations` field, along with their fixation time.  Further, "keys" to these
  mutations are removed from every gamete in the population.  This removal means that *selected mutations are no longer
  involved in genetic value calculations*. See :ref:`model_params` to point you to the docs on parameterizing
  simulations.

The "pruning" of fixations always results in faster simulations (assuming that you are simulating selected mutations).
Whether or not it is the right thing to do depends on the modeling scenario.  

Here ar examples of doing the wrong thing:

1. You are simulating additive effects on fitness, using :class:`fwdpy11.fitness.SlocusAdditive`.  If you prune
   fixations, you will be removing an *additive constant* from every diploid, meaning that relative fitnesses are
   affected.  We can illustrate this as follows, using a simple numerical example:

.. ipython:: python

    import numpy as np

    np.random.seed(42)
    # Pretend that these are the genetic values of 
    # ten diploids, only considering segregating variation:
    genetic_values = np.random.random_sample(10)
    # Pretend that this is the sum of effect sizes of fixations:
    fixed_effects = 0.1

    # Relative fitness, ignoring fixations:
    rw = genetic_values/np.sum(genetic_values)

    # Relative fitness, with fixations included
    rwf = (genetic_values+fixed_effects)/(np.sum(genetic_values+fixed_effects))
    
    # The values are different
    print(any(rw == rwf) is True)

    # But the rank orders are the same:
    rw_rank = np.argsort(rw)
    rwf_rank = np.argsort(rw)
    print(all(rw_rank == rwf_rank) is True)

2. You are simulating a trait under Gaussian stabilizing selection with respect to an optimum.  As mutations with effect
   sizes in the direction of the optimum increase in frequency, genotypes with those mutations increase in fitness.  If
   you prune selected fixations, each time a fixation occurs, the distance of every individual to the optimum is
   affected, and thus relative fitnesses may be affected.

General advice:

1. Prefer multiplicative fitness for regular pop-gen simulations.  Prune selected fixations in this case.
2. Do not prune selected fixations when simulating quantitative traits.

.. note::

    Neutral fixations are always recorded and removed from gametes.  This removal is for efficiency, so
    that we don't have to recombine fixed variants between gametes which, by definition, both have the fixed variants.

The future
-----------------------------------------------------

In the future, we hope to unify how fixations are handled.  Doing so will probably depend upon upstream changes in
fwdpp.
