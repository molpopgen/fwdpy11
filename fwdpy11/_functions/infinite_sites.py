from fwdpy11 import DiploidPopulation, GSLrng
from fwdpy11._fwdpy11 import _infinite_sites


def infinite_sites(rng: GSLrng, pop: DiploidPopulation, mutation_rate: float) -> int:
    return _infinite_sites(rng, pop, mutation_rate)
