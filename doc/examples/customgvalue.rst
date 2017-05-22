.. _customgvalues:

Custom stateless genetic value calculations
------------------------------------------------------------

Background:

* :ref:`definitions`

Further reading:

* For an example of a stateful model implemented in C++, see :ref:`stateful_fitness`.

fwdpy11 allows custom genetic value/fitness calculations to be implemented in either Python or C++.  The latter will
greatly outperform the former and we assume that Python implementations will be used for prototyping before writing
any C++ code.

Below, we will implement a pure additive fitness model with no dominance.  This will be equivalent to
:class:`fwdpy11.fitness.SlocusAdditive` with a scaling of 2 and a dominance of 1 for selected mutations.

We will evolve a population with our Python-based additive fitness model and compare timings to the built-in class,
which is implemented in C++.

.. ipython:: python
    :suppress:

    import fwdpy11
    import fwdpy11.fitness
    import fwdpy11.ezparams
    import fwdpy11.model_params
    import fwdpy11.python_genetic_values
    import fwdpy11.wright_fisher

This is our custom additive fitness function:

.. ipython:: python

    #This function suffices to calculate additive fitness
    def additive(dip,gametes,mutations):
        rv = sum([mutations[i].s for i in
            list(gametes[dip.first].smutations)+
            list(gametes[dip.second].smutations)])
        return max(0.0,1+rv)

We build a single-locus fitness calculator
by hooking our function up to an instance of 
:class:`fwdpy11.python_genetic_values.GeneticValue`:

.. ipython:: python

    custom_additive = fwdpy11.python_genetic_values.GeneticValue(additive)

Now, let's set up a population and get some model parameters 
established:

.. ipython:: python

    pop = fwdpy11.SlocusPop(1000)
    pdict = fwdpy11.ezparams.mslike(pop,
        dfe=fwdpy11.ExpS(0,1,1,-0.05),
        pneutral=0.95,simlen=10)
    pdict['gvalue'] = custom_additive
    params = fwdpy11.model_params.SlocusParams(**pdict)

    rng = fwdpy11.GSLrng(42)

Get mean run time using our custom additive model:

.. ipython:: python

    %timeit fwdpy11.wright_fisher.evolve(rng,pop,params)

OK, let's compare to the built-in additive fitness class.  We need
to reset our parameters, rng, etc., first:

.. ipython:: python

    pdict['gvalue'] = fwdpy11.fitness.SlocusAdditive(2.0)
    pop = fwdpy11.SlocusPop(1000)
    assert(pop.generation == 0)
    rng = fwdpy11.GSLrng(42)
    params = fwdpy11.model_params.SlocusParams(**pdict)

The C++ version is **much** faster:

.. ipython:: python

    %timeit fwdpy11.wright_fisher.evolve(rng,pop,params)

The reason for the massive speed difference has little to do with how we implemented our additive function. Rather, it
is due to a constant back and forth between C++ and Python. See `here
<https://pybind11.readthedocs.io/en/stable/advanced/cast/functional.html>`_ for details.  We are investigating ways to
eliminate the constant back and forth.  However, any changes will likely require refactoring some of the C++ code underlying this package, and so we're proceeding with caution.

.. todo::
    Provide example of a stateful model implemented in Python

.. _customgvaluecpp:

Implementing stateless genetic values in C++
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

In order to get our speed back for custom stateless genetic value calcuations, we need C++.  fwdpy11 provides a set of macros to greatly simplify implementing Python classes based on custom functions. *These will run at native C++ speed!*

Here is a complete example of a strictly additive effects on fitness model:

.. literalinclude:: custom_additive.cpp
    :language: cpp
    :lines: 19-

.. note::
    The defintion of the "additive" function looks very similar to our pur Python version above.  That is by design: fwdpy11 exposes the C++ types to Python directly, so that interacting with them involves the same concepts.

Let's use it:

.. ipython:: python

    import cppimport
    ca = cppimport.imp("doc.examples.custom_additive")
    pdict = fwdpy11.ezparams.mslike(pop,
        dfe=fwdpy11.ExpS(0,1,1,-0.05),
        pneutral=0.95,simlen=10)
    pdict['gvalue'] = ca.examples.custom_additive.additive()
    %timeit fwdpy11.wright_fisher.evolve(rng,pop,params)

That's much better--we are back to very fast fitness calculations!

.. note::
    In the above code block, the assignment of `pdict['gvalue']` appears to have an odd set of module names associated with it.  That is due
    to where the C++ source file is located for this example (in order that the manual build properly).  Typically, the .cpp file and the .py
    file are in the same directory, and the assignment would read `ca.additive()`.  See the unit test files tests/tests_custom_stateless_fitness.py
    for a more typical example.

Distinguishing Aa from aa genotypes
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The pure additive effects on fitness model is trivial to implement because it just requires summing over effect sizes.  In other words, the `h` field of a mutation is irrelevant.  In this section, I show an example of a fitness scheme where we must distinguish `Aa` from `aa` genotypes.  We will implement a model where the relative fitness of the three genotypes are :math:`1`, :math:`1+t`, and :math:`1+s`, respectively.  In order to do this, we will use the `h` field from a mutation for :math:`t` and the `s` field for :math:`s`.  For our example, the fitness of a diploid will be multiplicative over sites.

Here is our C++ code:

.. literalinclude:: custom_stateless_genotype.cpp
    :language: cpp
    :lines: 19-

How do we use the above? We are hijacking fwdpy11 machinery for our purposes.  In fwdpy11, new mutations have a fixed `h` and `s` can come from a distribution.  Thus, in our new model, `t` is fixed and `s` is either fixed or sampled from a distribution. Thus, when creating "regions" (:ref:`regions`) , you can get a variety of complex evolutionary scenarios happening.
