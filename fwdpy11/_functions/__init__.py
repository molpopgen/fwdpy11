from .data_matrix_from_tables import (data_matrix_from_tables,  # NOQA
                                      make_data_matrix)
from .import_demes import demography_from_demes  # NOQa
from .simplify_tables import simplify, simplify_tables  # NOQA

from fwdpy11._fwdpy11 import _infinite_sites
from fwdpy11 import GSLrng, DiploidPopulation


def infinite_sites(rng: GSLrng, pop: DiploidPopulation, mu: float) -> int:
    return _infinite_sites(rng, pop, mu)


def _validate_regions(regions, sequence_length):
    for r in regions:
        try:
            if r.beg >= sequence_length or r.end > sequence_length:
                raise ValueError(
                    f"Region {r} extends beyond the"
                    + f" sequence length {sequence_length}")
        except AttributeError:
            try:
                for des in r.des:
                    if des.beg >= sequence_length or des.end > sequence_length:
                        raise ValueError(
                            f"Region {des} in {r} extends beyond" +
                            f" the sequence length {sequence_length}")
            except TypeError:  # r.des not Iterable
                if r.des.beg >= sequence_length or r.des.end > sequence_length:
                    raise ValueError(
                        f"Region {des} in {r} extends beyond the" +
                        f" sequence length {sequence_length}")
