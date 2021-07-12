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
import demes
import fwdpy11
import numpy as np
import pytest


@pytest.fixture(scope="function")
def numpy_generator(request):
    try:
        seed = request.param.get("seed", None)
        return np.random.Generator(np.random.MT19937(seed=seed))
    except AttributeError as a:  # NOQA
        raise a

@pytest.fixture(scope="function")
def rng(request):
    try:
        seed = request.param.get("seed", None)
        return fwdpy11.GSLrng(seed)
    except AttributeError:  # NOQA
        return fwdpy11.GSLrng(42)


@pytest.fixture(scope="function")
def pop(request):
    popsize = request.param["N"]
    genome_length = request.param["genome_length"]
    return fwdpy11.DiploidPopulation(popsize, genome_length)


@pytest.fixture(scope="function")
def mslike_pop(request):
    try:
        N = request.param.get("N", None)
        return fwdpy11.DiploidPopulation(N, 1.0)
    except AttributeError:  # NOQA
        return fwdpy11.DiploidPopulation(1000, 1.0)


@pytest.fixture
def gutenkunst():
    yaml = """
description: The Gutenkunst et al. (2009) OOA model.
doi:
- https://doi.org/10.1371/journal.pgen.1000695
time_units: years
generation_time: 25

demes:
- name: ancestral
  description: Equilibrium/root population
  epochs:
  - {end_time: 220e3, start_size: 7300}
- name: AMH
  description: Anatomically modern humans
  ancestors: [ancestral]
  epochs:
  - {end_time: 140e3, start_size: 12300}
- name: OOA
  description: Bottleneck out-of-Africa population
  ancestors: [AMH]
  epochs:
  - {end_time: 21.2e3, start_size: 2100}
- name: YRI
  description: Yoruba in Ibadan, Nigeria
  ancestors: [AMH]
  epochs:
  - start_size: 12300
- name: CEU
  description: Utah Residents (CEPH) with Northern and Western European Ancestry
  ancestors: [OOA]
  epochs:
  - {start_size: 1000, end_size: 29725}
- name: CHB
  description: Han Chinese in Beijing, China
  ancestors: [OOA]
  epochs:
  - {start_size: 510, end_size: 54090}

migrations:
- {demes: [YRI, OOA], rate: 25e-5}
- {demes: [YRI, CEU], rate: 3e-5}
- {demes: [YRI, CHB], rate: 1.9e-5}
- {demes: [CEU, CHB], rate: 9.6e-5}
"""
    return demes.loads(yaml)


@pytest.fixture
def small_split_model():
    yaml = """
description: Two deme model.
time_units: generations
demes:
- name: ancestral
  description: ancestral deme
  epochs:
  - {start_size: 100}
- name: deme1
  description: child deme
  start_time: 100
  epochs:
  - {start_size: 100}
  ancestors: [ancestral]
"""
    return demes.loads(yaml)


@pytest.fixture
def two_deme_split_with_ancestral_size_change():
    yaml = """
description: Two deme model with migration and size changes.
time_units: generations
demes:
- name: ancestral
  description: ancestral deme, two epochs
  epochs:
  - {end_time: 20, start_size: 100}
  - {end_time: 10, start_size: 200}
- name: deme1
  description: child 1
  epochs:
  - {start_size: 250, end_size: 500, end_time: 0}
  ancestors: [ancestral]
- name: deme2
  description: child 2
  epochs:
  - {start_size: 50, end_size: 200, end_time: 0}
  ancestors: [ancestral]
migrations:
- {demes: [deme1, deme2], rate: 1e-3}
"""
    return demes.loads(yaml)


@pytest.fixture
def start_demes_model_with_two_pops():
    yaml = """
description: after slicing original graph at time 50 (below slice)
time_units: generations
demes:
- name: deme1
  epochs:
  - {start_size: 100, end_time: 50}
  - {start_size: 100, end_size: 200, end_time: 0}
- name: deme2
  epochs:
  - {start_size: 75, end_time: 50}
  - {start_size: 75, end_time: 20}
  - {start_size: 300, end_time: 0}
migrations:
- demes: [deme1, deme2]
  rate: 1e-2
  start_time: 50
"""
    return demes.loads(yaml)
