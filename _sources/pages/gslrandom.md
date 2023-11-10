(gslrandom)=

# Random number distributions

The following functions map directly to C functions
from the GNU Scientific Library.  All functions
take positional arguments without names.

```{eval-rst}
.. py:function:: fwdpy11.gsl_ran_gaussian_ziggurat

    :param rng: Random number generator
    :type rng: fwdpy11.GSLrng
    :param sd: The standard deviation
    :type sd: float
    :returns: A gaussian deviate with mean zero and standard
              deviation ``sd``
    :rtype: float

```

```{eval-rst}
.. py:function:: fwdpy11.gsl_rng_uniform

    :param rng: Random number generator
    :type rng: fwdpy11.GSLrng
    :returns: Uniform deviate from the range ``[0, 1)``
    :rtype: float

```

```{eval-rst}
.. py:function:: fwdpy11.gsl_ran_flat

    :param rng: Random number generator
    :type rng: fwdpy11.GSLrng
    :param a: Lower bound of range, inclusive.
    :type a: float
    :param b: Upper bound of range, exclusive.
    :type b: float
    :returns: Uniform deviate from the range ``[a, b)``.
    :rtype: float

```

```{eval-rst}
.. py:function:: fwdpy11.gsl_ran_poisson

    :param rng: Random number generator
    :type rng: fwdpy11.GSLrng
    :param mean: The mean of the distribution
    :type mean: float
    :returns: The number of successes
    :rtype: int

```

```{eval-rst}
.. py:function:: fwdpy11.gsl_ran_geometric

    :param rng: Random number generator
    :type rng: fwdpy11.GSLrng
    :param p: The probability of success
    :type p: float
    :returns: The waiting time until the next success
    :rtype: int

```


