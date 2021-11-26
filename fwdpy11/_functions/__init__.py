from .data_matrix_from_tables import (data_matrix_from_tables,  # NOQA
                                      make_data_matrix)
from .import_demes import demography_from_demes  # NOQa
from .simplify_tables import simplify, simplify_tables  # NOQA

from fwdpy11._fwdpy11 import _infinite_sites
from fwdpy11 import GSLrng, DiploidPopulation

def infinite_sites(rng: GSLrng, pop: DiploidPopulation, mu: float) -> int:
    return _infinite_sites(rng, pop, mu)
