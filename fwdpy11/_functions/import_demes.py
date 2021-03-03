from typing import Dict

import demes

from ..discrete_demography import (MassMigration, MigrationMatrix, SetDemeSize,
                                   SetExponentialGrowth, SetMigrationRates,
                                   SetSelfingRate, copy_individuals,
                                   move_individuals)


# TODO: need type hints for dg
def demography_from_demes(dg) -> Dict:
    demography: Dict = {
        "mass_migrations": [],
        "set_growth_rates": [],
        "set_deme_sizes": [],
        "set_selfing_rates": [],
        "migmatrix": None,
        "set_migration_rates": [],
    }
    return demography
