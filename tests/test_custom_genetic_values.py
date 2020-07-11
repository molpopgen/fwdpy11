#
# Copyright (C) 2017-2020 Kevin Thornton <krthornt@uci.edu>
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

import unittest

import numpy as np

import fwdpy11


class TestCustomGeneticValueisTrait(unittest.TestCase):
    def test_run(self):

        import pygss

        GSS = pygss.PyGSS(opt=0.0, VS=1.0)

        pdict = {
            "nregions": [],
            "sregions": [
                fwdpy11.mvDES(
                    fwdpy11.MultivariateGaussianEffects(
                        0, 1, 1, h=0.25, cov_matrix=np.identity(2)
                    ),
                    np.zeros(2),
                )
            ],
            "recregions": [fwdpy11.PoissonInterval(0, 1, 0.5)],
            "rates": (0, 1e-2, None),
            "demography": fwdpy11.DiscreteDemography(
                migmatrix=np.array([0.9, 0.1, 0.1, 0.9]).reshape((2, 2))
            ),
            "simlen": 100,
            "gvalue": fwdpy11.Additive(ndemes=2, scaling=2, gvalue_to_fitness=GSS,),
        }

        params = fwdpy11.ModelParams(**pdict)
        pop = fwdpy11.DiploidPopulation([1000, 1000], 1.0)

        rng = fwdpy11.GSLrng(1010)
        fwdpy11.evolvets(rng, pop, params, 100, suppress_table_indexing=True)

        self.assertEqual(pop.generation, params.simlen)
        md = np.array(pop.diploid_metadata)

        pdict["gvalue"] = fwdpy11.Additive(
            # ndemes=2, scaling=2, gvalue_to_fitness=
            ndemes=2,
            scaling=2,
            gvalue_to_fitness=fwdpy11.GSS(fwdpy11.Optimum(optimum=0.0, VS=1.0)),
        )
        params = fwdpy11.ModelParams(**pdict)
        pop2 = fwdpy11.DiploidPopulation([1000, 1000], 1.0)
        rng = fwdpy11.GSLrng(1010)
        fwdpy11.evolvets(rng, pop2, params, 100, suppress_table_indexing=True)
        md2 = np.array(pop.diploid_metadata)
        self.assertAlmostEqual(md["g"].mean(), md2["g"].mean())

    def test_run_stateful(self):

        import pygss

        GSS = pygss.PyGSSRandomOptimum(opt=0.0, VS=1.0)

        pdict = {
            "nregions": [],
            "sregions": [
                fwdpy11.mvDES(
                    fwdpy11.MultivariateGaussianEffects(
                        0, 1, 1, h=0.25, cov_matrix=np.identity(2)
                    ),
                    np.zeros(2),
                )
            ],
            "recregions": [fwdpy11.PoissonInterval(0, 1, 0.5)],
            "rates": (0, 1e-2, None),
            "demography": fwdpy11.DiscreteDemography(
                migmatrix=np.array([0.9, 0.1, 0.1, 0.9]).reshape((2, 2))
            ),
            "simlen": 5,
            "gvalue": fwdpy11.Additive(ndemes=2, scaling=2, gvalue_to_fitness=GSS,),
        }

        params = fwdpy11.ModelParams(**pdict)
        pop = fwdpy11.DiploidPopulation([1000, 1000], 1.0)

        rng = fwdpy11.GSLrng(1010)
        fwdpy11.evolvets(rng, pop, params, 100, suppress_table_indexing=True)
        self.assertEqual(
            len(params.gvalue.gvalue_to_fitness.optima), pop.generation + 1
        )


class TestCustomGeneticValueNoise(unittest.TestCase):
    def test_run(self):

        import pynoise

        GSS = fwdpy11.GSS(optimum=0.0, VS=1.0)
        Noise = pynoise.PyNoise()

        pdict = {
            "nregions": [],
            "sregions": [
                fwdpy11.mvDES(
                    fwdpy11.MultivariateGaussianEffects(
                        0, 1, 1, h=0.25, cov_matrix=np.identity(2)
                    ),
                    np.zeros(2),
                )
            ],
            "recregions": [fwdpy11.PoissonInterval(0, 1, 0.5)],
            "rates": (0, 1e-2, None),
            "demography": fwdpy11.DiscreteDemography(
                migmatrix=np.array([0.9, 0.1, 0.1, 0.9]).reshape((2, 2))
            ),
            "simlen": 100,
            "gvalue": fwdpy11.Additive(
                ndemes=2, scaling=2, gvalue_to_fitness=GSS, noise=Noise
            ),
        }

        params = fwdpy11.ModelParams(**pdict)
        pop = fwdpy11.DiploidPopulation([1000, 1000], 1.0)

        rng = fwdpy11.GSLrng(1010)
        fwdpy11.evolvets(rng, pop, params, 100, suppress_table_indexing=True)

        self.assertEqual(pop.generation, params.simlen)
        md = np.array(pop.diploid_metadata)
        self.assertTrue(len(np.where(md['e'] > 0)[0]) > 0)


if __name__ == "__main__":
    unittest.main()
