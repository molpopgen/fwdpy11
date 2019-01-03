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

import sys

if sys.version_info[0] < 3:
    raise ValueError("Python3 required!")

from fwdpy11._init import * # NOQA
from fwdpy11._version import __version__ # NOQA
from ._fwdpy11 import * # NOQA

from .fwdpp_types import *
from ._opaque_gametes import *
from ._opaque_mutations import *
from ._opaque_diploids import *
from ._Population import VecUint32
from .fwdpy11_types import DiploidGenotype
from .fwdpy11_types import DiploidMetadata
from .fwdpy11_types import Mutation
# from ._regions import *
from ._dev import *
from ._gslrng import GSLrng
from ._Population import Population
from .SlocusPop import SlocusPop
from .MlocusPop import MlocusPop

