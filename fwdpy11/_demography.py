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
    return np.exp((np.log(Nt) - np.log(N0))/time)
