#
# Copyright (C) 2023 Kevin Thornton <krthornt@uci.edu>
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
import pytest

import fwdpy11


def _validate_expected(regions, numbers):
    assert len([i for i in regions if isinstance(
        i, fwdpy11.PoissonCrossoverGenerator)]) == numbers[0]

    assert len([i for i in regions if isinstance(
        i, fwdpy11.NonPoissonCrossoverGenerator)]) == numbers[1]


def roundtrip(regions):
    yaml = """
    time_units: generations
    demes:
    - name: A
      epochs:
        - {end_time: 0, start_size: 10}
    """
    graph = demes.loads(yaml)
    demog = fwdpy11.discrete_demography.from_demes(graph, burnin=1)
    pdict = {"recregions": regions,
             "nregions": [],
             "sregions": [],
             "rates": (0, 0, None),
             "gvalue": fwdpy11.Multiplicative(2.),
             "simlen": 10,
             "demography": demog,
             }
    params = fwdpy11.ModelParams(**pdict)
    pop = fwdpy11.DiploidPopulation(10, 2)
    rng = fwdpy11.GSLrng(123)
    fwdpy11.evolvets(rng, pop, params, 100)


# The setup is ([regions], (number poisson, number non-poisson))
# The latter is used to mimic the internal idiom that we
# use to set up the genetic map to send to C++
@pytest.mark.parametrize("data", [
    ([fwdpy11.PoissonInterval(0, 1, 1e-3)], (1, 0)),
    ([fwdpy11.PoissonPoint(1, 1e-3)], (1, 0)),
    ([fwdpy11.PoissonInterval(0, 1, 1e-3),
      fwdpy11.PoissonPoint(1, 1e-3)], (2, 0)),
    ([fwdpy11.BinomialInterval(0, 1, 1e-3)], (0, 1)),
    ([fwdpy11.BinomialPoint(1, 1e-3)], (0, 1)),
    ([fwdpy11.BinomialInterval(0, 1, 1e-3),
      fwdpy11.BinomialPoint(1, 1e-3)], (0, 2)),
    ([fwdpy11.PoissonInterval(0, 0.1, 1e-3),
     fwdpy11.PoissonPoint(0.1, 1e-3),
     fwdpy11.BinomialInterval(0.1, 0.2, 1e-3),
     fwdpy11.BinomialPoint(0.2, 1e-3)], (2, 2)),
    ([fwdpy11.BinomialIntervalMap(0.5, [fwdpy11.Region(0, 1, 1e-5),
                                        fwdpy11.Region(1, 2, 2)])], (0, 1))
])
def test_model_params_with_recregions(data):
    regions, expected = data
    _validate_expected(regions, expected)
    roundtrip(regions)
