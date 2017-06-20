import math


def exponential_size_change(Nstart, Nstop, time):
    """
    Generate a list of population sizes
    according to exponential size_change model

    :param Nstart: population size at onset of size change
    :param Nstop: Population size to reach at end of size change
    :param time: Time (in generations) to get from Nstart to Nstop

    :return: A list of integers representing population size over time.

    .. versionadded:: 0.1.1
    """
    if time < 1:
        raise RuntimeError("time must be >= 1")
    if Nstart < 1 or Nstop < 1:
        raise RuntimeError("Nstart and Nstop must both be >= 1")
    G = math.exp((math.log(Nstop) - math.log(Nstart))/time)
    rv = []
    for i in range(time):
        rv.append(round(Nstart*pow(G, i+1)))
    return rv
