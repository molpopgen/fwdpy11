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


def build_model(gvalue_to_fitness, noise, simlen, sizes):
    yaml = f"""
time_units: generations
demes:
 - name: A
   epochs:
    - start_size: {sizes[0]}
 - name: B
   epochs:
    - start_size: {sizes[1]}
migrations:
 - demes: [A, B]
   rate: 0.1
"""
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
        "demography": fwdpy11.ForwardDemesGraph.from_demes(
            yaml, burnin=simlen, burnin_is_exact=True
        ),
        "simlen": simlen,
        "gvalue": fwdpy11.Additive(
            ndemes=2,
            scaling=2,
            gvalue_to_fitness=gvalue_to_fitness,
            noise=noise,
        ),
    }
    return pdict


def build_single_trait_model(gvalue, simlen):
    pdict = {
        "nregions": [],
        "sregions": [fwdpy11.GaussianS(0.0, 1.0, 1.0, 0.1)],
        "recregions": [fwdpy11.PoissonInterval(0, 1, 0.5)],
        "rates": (0, 1e-3, None),
        "simlen": simlen,
        "gvalue": gvalue,
    }
    return pdict


class TestCustomGeneticValueisTrait(unittest.TestCase):
    def test_run(self):
        import pygss

        GSS = pygss.PyGSS(opt=0.0, VS=1.0)
        pdict = build_model(GSS, None, 100, [1000, 1000])
        pdict["prune_selected"] = False

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
            gvalue_to_fitness=fwdpy11.GaussianStabilizingSelection.single_trait(
                [fwdpy11.Optimum(optimum=0.0, VS=1.0)]
            ),
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
        pdict = build_model(GSS, None, 5, [1000, 1000])
        pdict["prune_selected"] = False
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

        GSS = fwdpy11.fwdpy11.GaussianStabilizingSelection.single_trait(
            [fwdpy11.Optimum(optimum=0.0, VS=1.0, when=0)]
        )
        Noise = pynoise.PyNoise()
        pdict = build_model(GSS, Noise, 5, [1000, 1000])
        pdict["prune_selected"] = False
        params = fwdpy11.ModelParams(**pdict)
        pop = fwdpy11.DiploidPopulation([1000, 1000], 1.0)

        rng = fwdpy11.GSLrng(1010)
        fwdpy11.evolvets(rng, pop, params, 100, suppress_table_indexing=True)

        self.assertEqual(pop.generation, params.simlen)
        md = np.array(pop.diploid_metadata)
        self.assertTrue(len(np.where(md["e"] > 0)[0]) > 0)


class TestCustomPyGeneticValue(unittest.TestCase):
    @classmethod
    def setUp(self):
        import pygss
        import pyadditive
        import pynoise

        PyA = pyadditive.PyAdditive(
            pygss.PyGSS(opt=0.0, VS=1.0), noise=pynoise.PyNoise()
        )

        self.pdict_PyA = build_single_trait_model(PyA, 100)

    def test_run(self):
        pop_PyA = fwdpy11.DiploidPopulation(1000, 1.0)
        self.pdict_PyA["demography"] = fwdpy11.ForwardDemesGraph.tubes(
            pop_PyA.deme_sizes()[1],
            burnin=self.pdict_PyA["simlen"],
            burnin_is_exact=True,
        )
        self.pdict_PyA["prune_selected"] = False
        mparams_PyA = fwdpy11.ModelParams(**self.pdict_PyA)
        rng_PyA = fwdpy11.GSLrng(42 * 666)

        fwdpy11.evolvets(rng_PyA, pop_PyA, mparams_PyA, 100)
        if len(pop_PyA.tables.mutations) == 0:
            self.fail("simulation ended with no mutations")
        md = np.array(pop_PyA.diploid_metadata, copy=False)
        self.assertTrue(md["g"].var() > 0.0)
        self.assertTrue(md["e"].var() > 0.0)


class TestCustomPyGeneticValueOverloadGeneticValueToFitness(unittest.TestCase):
    def build_model(self, gvalue, simlen):
        pdict = {
            "nregions": [],
            "sregions": [fwdpy11.GaussianS(beg=0.0, end=1.0, weight=1.0, sd=0.1)],
            "recregions": [fwdpy11.PoissonInterval(0, 1, 0.5)],
            "rates": (0, 1e-2, None),
            "demography": None,
            "simlen": simlen,
            "gvalue": gvalue,
            "prune_selected": False,
        }
        return pdict

    def test_run(self):
        import pyadditivegss

        pygv = pyadditivegss.PyAdditiveGSS(opt=0.0, VS=1.0)
        pdict = self.build_model(pygv, 100)
        pop = fwdpy11.DiploidPopulation(1000, 1.0)
        pdict["demography"] = fwdpy11.ForwardDemesGraph.tubes(
            pop.deme_sizes()[1], burnin=pdict["simlen"], burnin_is_exact=True
        )

        params = fwdpy11.ModelParams(**pdict)

        rng = fwdpy11.GSLrng(1010)
        fwdpy11.evolvets(rng, pop, params, 100, suppress_table_indexing=True)

        self.assertEqual(pop.generation, params.simlen)
        md = np.array(pop.diploid_metadata)

        pdict["gvalue"] = fwdpy11.Additive(
            2.0,
            fwdpy11.GaussianStabilizingSelection.single_trait(
                [fwdpy11.Optimum(optimum=0.0, VS=1.0, when=0)]
            ),
        )
        params = fwdpy11.ModelParams(**pdict)
        pop2 = fwdpy11.DiploidPopulation(1000, 1.0)
        rng = fwdpy11.GSLrng(1010)
        fwdpy11.evolvets(rng, pop2, params, 100, suppress_table_indexing=True)
        md2 = np.array(pop.diploid_metadata)
        self.assertAlmostEqual(md["g"].mean(), md2["g"].mean())


if __name__ == "__main__":
    unittest.main()
