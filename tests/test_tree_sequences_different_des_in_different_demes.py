#
# Copyright (C) 2020 Kevin Thornton <krthornt@uci.edu>
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

# This file contains round-trip tests of correctness using
# fwdpy11.mvDES to specify different DES in different demes.
# Low-level tests of fwdpy11.mvDES itself are in test_regions.py

import unittest
import pytest

import demes
import fwdpy11
import numpy as np


def gvalue_multiplicative(pop, ind, scaling):
    g = 1.0
    keys = [k for k in pop.haploid_genomes[pop.diploids[ind].first].smutations]
    keys.extend(
        [k for k in pop.haploid_genomes[pop.diploids[ind].second].smutations])
    keys = np.array(keys, dtype=np.uint32)
    key_counts = np.unique(keys, return_counts=True)

    ind_md = pop.diploid_metadata[ind]
    for i, j in zip(key_counts[0], key_counts[1]):
        if j == 1:
            g *= (
                1.0
                + pop.mutations[i].esizes[ind_md.deme]
                * pop.mutations[i].heffects[ind_md.deme]
            )
        elif j == 2:
            g *= 1.0 + scaling * pop.mutations[i].esizes[ind_md.deme]

    return g


class TestMultiplicativeWithExpSNoMigration(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.pop = fwdpy11.DiploidPopulation([50, 50], 1)
        mvDES = fwdpy11.mvDES(
            [fwdpy11.ExpS(0, 1, 1, 0.1), fwdpy11.ExpS(0, 1, 1, -0.1)],
            np.zeros(2),
            np.identity(2),
        )

        self.pdict = {
            "nregions": [],
            "sregions": [mvDES],
            "recregions": [fwdpy11.PoissonInterval(0, 1, 1e-2)],
            "rates": (0, 1e-2, None),
            "gvalue": fwdpy11.Multiplicative(2.0, ndemes=2),
            "demography": fwdpy11.ForwardDemesGraph.tubes(self.pop.deme_sizes()[1], 1),
            "simlen": 100,
            "prune_selected": True,
        }
        self.params = fwdpy11.ModelParams(**self.pdict)
        self.rng = fwdpy11.GSLrng(918273)
        fwdpy11.evolvets(self.rng, self.pop, self.params, 100)
        n0 = 0
        n1 = 0
        for i, md in enumerate(self.pop.diploid_metadata):
            ng0 = len(
                self.pop.haploid_genomes[self.pop.diploids[i].first].smutations)
            ng1 = len(
                self.pop.haploid_genomes[self.pop.diploids[i].second].smutations)
            if ng0 + ng1 > 0:
                if md.deme == 0:
                    n0 += 1
                elif md.deme == 1:
                    n1 += 1
        if n0 == 0 or n1 == 0:
            self.assertFail(
                "we don't have individuals with mutations in each deme")

    def test_it(self):
        gv = np.zeros(self.pop.N)
        for i in range(self.pop.N):
            gv[i] = gvalue_multiplicative(self.pop, i, 2.0)
        md = np.array(self.pop.diploid_metadata, copy=False)
        self.assertTrue(np.allclose(gv, md["g"]))


class TestGaussianStabilizingSelection(unittest.TestCase):
    """
    Quick sim w/high migration rate
    to put the same mutation in both demes.
    """

    @ classmethod
    def setUpClass(self):
        yaml = """
        time_units: generations
        demes:
         - name: A
           epochs:
            - start_size: 100
         - name: B
           epochs:
            - start_size: 100
        migrations:
         - demes: [A, B]
           rate: 0.1
        """
        demography = fwdpy11.ForwardDemesGraph.from_demes(yaml, 1)
        pdict = {
            "nregions": [],
            "sregions": [
                fwdpy11.mvDES(
                    fwdpy11.MultivariateGaussianEffects(
                        0, 1, 1, h=1, cov_matrix=np.identity(2)
                    ),
                    np.zeros(2),
                )
            ],
            "recregions": [],
            "rates": (0, 1e-3, None),
            "demography": demography,
            "simlen": 100,
            "gvalue": fwdpy11.Additive(
                ndemes=2, scaling=2, gvalue_to_fitness=fwdpy11.GSS(optimum=0.0, VS=1.0)
            ),
        }

        self.params = fwdpy11.ModelParams(**pdict)
        self.pop = fwdpy11.DiploidPopulation([100, 100], 1.0)
        self.rng = fwdpy11.GSLrng(1010)
        fwdpy11.evolvets(self.rng, self.pop, self.params, 10)

        tv = fwdpy11.TreeIterator(
            self.pop.tables, self.pop.alive_nodes, update_samples=True
        )
        nt = np.array(self.pop.tables.nodes, copy=False)
        demes_with_muts = np.zeros(2)
        for t in tv:
            for m in t.mutations():
                demes = np.unique(nt["deme"][t.samples_below(m.node)])
                for d in demes:
                    demes_with_muts[d] += 1
        assert np.all(demes_with_muts >
                      0), "test requires mutations in both demes"

    def test_genetic_values(self):
        for m in self.pop.diploid_metadata:
            g = 0.0
            for i in [
                self.pop.diploids[m.label].first,
                self.pop.diploids[m.label].second,
            ]:
                for k in self.pop.haploid_genomes[i].smutations:
                    g += self.pop.mutations[k].esizes[m.deme]
            self.assertAlmostEqual(m.g, g)


class TestMultivariateLogNormalS(unittest.TestCase):
    """
    Quick sim w/high migration rate
    to put the same mutation in both demes.
    """

    @ classmethod
    def setUpClass(self):
        yaml = """
        time_units: generations
        demes:
         - name: A
           epochs:
            - start_size: 100
         - name: B
           epochs:
            - start_size: 100
        migrations:
         - demes: [A, B]
           rate: 0.1
        """
        demography = fwdpy11.ForwardDemesGraph.from_demes(yaml, 1)
        ln = fwdpy11.LogNormalS.mv(0, 1, 1, scaling=-200)
        pdict = {
            "nregions": [],
            "sregions": [fwdpy11.mvDES(ln, np.zeros(2), np.identity(2))],
            "recregions": [],
            "rates": (0, 1e-3, None),
            "demography": demography,
            "simlen": 100,
            "gvalue": fwdpy11.Multiplicative(ndemes=2, scaling=2),
        }

        self.params = fwdpy11.ModelParams(**pdict)
        self.pop = fwdpy11.DiploidPopulation([100, 100], 1.0)
        self.rng = fwdpy11.GSLrng(201012)
        fwdpy11.evolvets(self.rng, self.pop, self.params, 10)

        tv = fwdpy11.TreeIterator(
            self.pop.tables, self.pop.alive_nodes, update_samples=True
        )
        nt = np.array(self.pop.tables.nodes, copy=False)
        demes_with_muts = np.zeros(2)
        for t in tv:
            for m in t.mutations():
                demes = np.unique(nt["deme"][t.samples_below(m.node)])
                for d in demes:
                    demes_with_muts[d] += 1
        assert np.all(demes_with_muts >
                      0), "test requires mutations in both demes"

    def test_genetic_values(self):
        for i, m in enumerate(self.pop.diploid_metadata):
            g = gvalue_multiplicative(self.pop, i, 2.0)
            self.assertAlmostEqual(m.g, g)


@ pytest.mark.parametrize("mvDES",
                          [fwdpy11.mvDES(fwdpy11.MultivariateGaussianEffects(0, 1, 1, np.identity(2)), np.zeros(2)),
                           fwdpy11.mvDES(
                              [fwdpy11.ConstantS(0, 1, 1, 0.1),
                               fwdpy11.ConstantS(0, 1, 1, -0.1)],
                              np.zeros(2),
                              np.identity(2),
                          )])
@ pytest.mark.parametrize("gvalue", [fwdpy11.Multiplicative(2., ndemes=2), fwdpy11.Additive(2., ndemes=2)])
def test_invalid_model(mvDES, gvalue):
    demog = """
description: trigger exception
time_units: generations
demes:
  - name: ancestor
    epochs:
      - start_size: 100
        end_time: 100
  - name: A
    ancestors: [ancestor]
    epochs:
      - start_size: 100
  - name: B
    ancestors: [ancestor]
    epochs:
      - start_size: 100
  - name: C
    ancestors: [ancestor]
    epochs:
      - start_size: 100
    """

    g = demes.loads(demog)
    model = fwdpy11.discrete_demography.from_demes(g, burnin=1)
    pdict = {
        "nregions": [],
        "sregions": [mvDES],
        "recregions": [fwdpy11.PoissonInterval(0, 1, 1e-2)],
        "rates": (0, 1, None),
        "gvalue": gvalue,
        "demography": model,
        "simlen": model.metadata["total_simulation_length"],
        "prune_selected": True,
    }
    params = fwdpy11.ModelParams(**pdict)
    rng = fwdpy11.GSLrng(918273)
    pop = fwdpy11.DiploidPopulation(100, 1.0)
    with pytest.raises(ValueError):
        fwdpy11.evolvets(rng, pop, params, 100)


if __name__ == "__main__":
    unittest.main()
