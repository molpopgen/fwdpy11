"""
tskit flags are unsigned 32 bit integers.

The first 16 bits are reserved for internal
use by tskit, so 16 and on are available
to us.

See here for details:
    https://tskit.readthedocs.io/en/latest/data-model.html?highlight=reserved#table-definitions

"""
INDIVIDUAL_IS_ALIVE = 2 ** 16
INDIVIDUAL_IS_PRESERVED = 2 ** 17
INDIVIDUAL_IS_FIRST_GENERATION = 2 ** 18
