def exponential_growth_rate(N0, Nt, time):
    """
    :param N0: Population size at onset of size change
    :param Nt: Population size to reach at end of size change
    :param time: Time (in generations) to get from N0 to Nt

    :return: A list of integers representing population size over time.

    .. versionadded:: 0.6.0
    """
    if time < 1:
        raise ValueError("time must be >= 1")
    if N0 < 1 or Nt < 1:
        raise ValueError("N0 and Nt must both be >= 1")
    import numpy as np

    return np.exp((np.log(Nt) - np.log(N0)) / time)


def migration_matrix_single_extant_deme(ndemes, focal_deme):
    """
    Sets up a migration matrix for a system of ndemes
    but where only one deme has a size > 0.  This
    is a convenience function for simulations where
    demes will appear later on and have migration amongst them.

    :param ndemes: The number of demes that will eventually exist.
    :param focal_deme: The index of the initial extant deme.
    :rtype: numpy.ndarray

    .. versionadded:: 0.6.0
    """
    import numpy as np

    m = np.zeros(ndemes * ndemes).reshape(ndemes, ndemes)
    m[focal_deme, focal_deme] = 1.0
    return m
