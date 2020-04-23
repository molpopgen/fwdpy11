#
# Copyright (C) 2017 Kevin Thornton <krthornt@uci.edu>
#
# This file is part of fwdpye1.
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

# This is part of fwdpy11's test suite.
# This file contains functions for quickly
# obtaining simulated populations


def quick_neutral_slocus(N=1000, simlen=100):
    from fwdpy11.ezparams import mslike
    from fwdpy11 import ModelParams
    from fwdpy11 import DiploidPopulation, GSLrng
    from fwdpy11 import Multiplicative
    from fwdpy11 import evolve_genomes

    pop = DiploidPopulation(N)
    params_dict = mslike(pop, simlen=simlen)
    params_dict["gvalue"] = Multiplicative(2.0)
    params = ModelParams(**params_dict)
    rng = GSLrng(42)
    evolve_genomes(rng, pop, params)
    return pop


def quick_nonneutral_slocus(N=1000, simlen=100, dfe=None):
    from fwdpy11.ezparams import mslike
    from fwdpy11 import ModelParams
    from fwdpy11 import DiploidPopulation, GSLrng
    from fwdpy11 import evolve_genomes
    from fwdpy11 import ExpS
    from fwdpy11 import Multiplicative

    pop = DiploidPopulation(N)
    if dfe is None:
        dfe = ExpS(0, 1, 1, -0.1)
    params_dict = mslike(pop, simlen=simlen, dfe=dfe, pneutral=0.95)
    params_dict["gvalue"] = Multiplicative(2.0)
    params = ModelParams(**params_dict)
    rng = GSLrng(42)
    evolve_genomes(rng, pop, params)
    return pop


def quick_slocus_qtrait_pop_params(N=1000, simlen=100):
    from fwdpy11 import DiploidPopulation
    from fwdpy11 import GSS
    from fwdpy11 import Additive
    from fwdpy11 import GaussianS, Region
    import numpy as np

    p = {
        "nregions": [],
        "sregions": [GaussianS(0, 1, 1, 0.25)],
        "recregions": [Region(0, 1, 1)],
        "rates": (0.0, 2e-3, 1e-3),
        "popsizes": np.array([N] * simlen, dtype=np.uint32),
        "simlen": simlen,
        "gvalue": Additive(2.0, GSS(VS=1.0, optimum=0.0)),
        "prune_selected": False,
    }
    pop = DiploidPopulation(N)
    return (pop, p)
