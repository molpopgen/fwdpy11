#
# Copyright (C) 2017 Kevin Thornton <krthornt@uci.edu>
#
# This file is part of fwdpy11.
#
# fwdpy11 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# fwdpy11 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with fwdpy11.  If not, see <http://www.gnu.org/licenses/>.
#

""" fwdpy11

    isort:skip_file
"""
import sys

from fwdpy11._version import version as __version__  # NOQA

from ._demography import *  # NOQA
from ._dev import *  # NOQA
from ._fwdpy11 import *  # NOQA
from ._monkeypatch import _data_matrix  # NOQA
from ._monkeypatch import _diploid_population  # NOQA
from ._monkeypatch import _table_collection  # NOQA
from .discrete_demography import (  # NOQA
    DiscreteDemography,
    MassMigration,
    move_individuals,
    copy_individuals,
    SetDemeSize,
    MigrationMatrix,
    SetExponentialGrowth,
    SetMigrationRates,
    SetSelfingRate,
)
from .regions import *  # NOQA
from .genetic_map_unit import (
    PoissonPoint,
    PoissonInterval,
    BinomialInterval,
    BinomialPoint,
    FixedCrossovers,
)
from .genetic_values import (  # NOQA
    PleiotropicOptima,
    Optimum,
    GSS,
    GSSmo,
    MultivariateGSS,
    MultivariateGSSmo,
    NoNoise,
    GaussianNoise,
    Additive,
    Multiplicative,
    GBR,
    StrictAdditiveMultivariateEffects,
)
from ._types.demography_debugger import DemographyDebugger  # NOQA
from ._types.model_params import ModelParams, MutationAndRecombinationRates  # NOQA
from ._evolve_genomes import *  # NOQA
from ._evolvets import *  # NOQA

if sys.version_info[0] < 3:
    raise ValueError("Python3 required!")


# NOTE: some operations that can be implemented efficiently
# in Python are supplied as monkey-patches to the pybind11 classes
_monkeypatch._diploid_population._patch_diploid_population(DiploidPopulation)  # NOQA
_monkeypatch._table_collection._patch_table_collection(TableCollection)  # NOQA
_monkeypatch._data_matrix._patch_data_matrix(DataMatrix)  # NOQA
