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
    from fwdpy11.model_params import ModelParams
    from fwdpy11 import SlocusPop, GSLrng
    from fwdpy11.genetic_values import SlocusMult
    from fwdpy11.wright_fisher import evolve
    pop = SlocusPop(N)
    params_dict = mslike(pop, simlen=simlen)
    params_dict['gvalue'] = SlocusMult(2.)
    params = ModelParams(**params_dict)
    rng = GSLrng(42)
    evolve(rng, pop, params)
    return pop


def quick_nonneutral_slocus(N=1000, simlen=100, dfe=None):
    from fwdpy11.ezparams import mslike
    from fwdpy11.model_params import ModelParams
    from fwdpy11 import SlocusPop, GSLrng
    from fwdpy11.wright_fisher import evolve
    from fwdpy11 import ExpS
    from fwdpy11.genetic_values import SlocusMult
    pop = SlocusPop(N)
    if dfe is None:
        dfe = ExpS(0, 1, 1, -0.1)
    params_dict = mslike(
        pop, simlen=simlen, dfe=dfe,
        pneutral=0.95)
    params_dict['gvalue'] = SlocusMult(2.0)
    params = ModelParams(**params_dict)
    rng = GSLrng(42)
    evolve(rng, pop, params)
    return pop


def quick_slocus_qtrait_pop_params(N=1000, simlen=100):
    from fwdpy11 import SlocusPop
    from fwdpy11.genetic_values import GSS
    from fwdpy11.genetic_values import SlocusAdditive
    from fwdpy11 import GaussianS, Region
    import numpy as np
    p = {'nregions': [],
         'sregions': [GaussianS(0, 1, 1, 0.25)],
         'recregions': [Region(0, 1, 1)],
         'rates': (0.0, 2e-3, 1e-3),
         'demography': np.array([N] * simlen, dtype=np.uint32),
         'gvalue': SlocusAdditive(2.0, GSS(VS=1.0, opt=0.0)),
         'prune_selected': False
         }
    pop = SlocusPop(N)
    return (pop, p)


def quick_mlocus_qtrait_pop_params(N=1000, simlen=100):
    from fwdpy11 import MlocusPop
    from fwdpy11.genetic_values import GSS
    from fwdpy11.genetic_values import MlocusAdditive
    from fwdpy11 import GaussianS, Region
    from fwdpy11.multilocus import binomial_rec
    import numpy as np

    theta, rho = 100., 100.
    mu = 1e-3
    sigmu = 0.25
    nloci = 5
    # Simulate a nloci locus system where each locus is 10 units long
    locus_boundaries = [(float(i + i * 10), float(i + i * 10 + 10))
                        for i in range(nloci)]
    nregions = [[Region(j[0], j[1], theta / (4. * float(N)), coupled=True)]
                for i, j in zip(range(nloci), locus_boundaries)]
    recregions = [[Region(j[0], j[1], rho / (4. * float(N)), coupled=True)]
                  for i, j in zip(range(nloci), locus_boundaries)]
    sregions = [[GaussianS(j[0] + 5., j[0] + 6., mu, sigmu, coupled=False)]
                for i, j in zip(range(nloci), locus_boundaries)]
    interlocus_rec = binomial_rec([0.5] * (nloci - 1))
    nlist = np.array([N] * simlen, dtype=np.uint32)
    mutrates_n = [theta / (4. * float(N))] * nloci
    mutrates_s = [mu] * nloci
    recrates = [rho / (4. * float(N))] * nloci
    param_dict = {'nregions': nregions,
                  'sregions': sregions,
                  'recregions': recregions,
                  'rates': (mutrates_n, mutrates_s, recrates),
                  'interlocus_rec': interlocus_rec,
                  'gvalue': MlocusAdditive(2.0, GSS(VS=1.0, opt=0.0)),
                  'demography': nlist,
                  'prune_selected': False
                  }
    pop = MlocusPop(N, locus_boundaries)
    return (pop, param_dict)


def quick_mlocus_qtrait(N=1000, simlen=100):
    from fwdpy11 import GSLrng
    from fwdpy11.wright_fisher import evolve
    from fwdpy11.model_params import ModelParams
    rng = GSLrng(42)
    pop, param_dict = quick_mlocus_qtrait_pop_params(N, simlen)
    params = ModelParams(**param_dict)
    evolve(rng, pop, params)
    return pop


# TODO: fix or remove
def quick_mlocus_qtrait_change_optimum(N=1000,
                                       simlen=100, prune_selected=False):
    """
    .. warning:: May result in long-running tests
    """
    from fwdpy11.model_params import ModelParams
    from fwdpy11.genetic_values import MlocusAdditive
    from fwdpy11.genetic_values import GSSmo
    from fwdpy11 import MlocusPop, GSLrng
    from fwdpy11.wright_fisher import evolve
    from fwdpy11 import GaussianS, Region
    from fwdpy11.multilocus import binomial_rec
    import numpy as np

    rng = GSLrng(42)
    theta, rho = 100., 100.
    mu = 1e-3
    sigmu = 0.25
    nloci = 5
    # Simulate a nloci locus system where each locus is 10 units long
    locus_boundaries = [(float(i + i * 10), float(i + i * 10 + 10))
                        for i in range(nloci)]
    nregions = [[Region(j[0], j[1], theta / (4. * float(N)), coupled=True)]
                for i, j in zip(range(nloci), locus_boundaries)]
    recregions = [[Region(j[0], j[1], rho / (4. * float(N)), coupled=True)]
                  for i, j in zip(range(nloci), locus_boundaries)]
    sregions = [[GaussianS(j[0] + 5., j[0] + 6., mu, sigmu, coupled=False)]
                for i, j in zip(range(nloci), locus_boundaries)]
    interlocus_rec = binomial_rec([0.5] * (nloci - 1))
    nlist = np.array([N] * simlen, dtype=np.uint32)
    param_dict = {'nregions': nregions,
                  'sregions': sregions,
                  'recregions': recregions,
                  'interlocus': interlocus_rec,
                  'mutrates_n': [theta / (4. * float(N))] * nloci,
                  'mutrates_s': [mu] * nloci,
                  'recrates': [rho / (4. * float(N))] * nloci,
                  'gvalue': MlocusAdditive(2.0, GSSmo([(0, 0, 1),
                                                       (simlen / 2, 1, 1)])),
                  'demography': nlist,
                  'prune_selected': prune_selected}
    params = ModelParams(**param_dict)
    pop = MlocusPop(N, locus_boundaries)
    evolve(rng, pop, params)
    return pop
