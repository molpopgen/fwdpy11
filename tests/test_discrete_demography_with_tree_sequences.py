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
import pickle
import unittest
from collections import namedtuple

import numpy as np

import fwdpy11


def validate_alive_node_metadata(pop):
    nodes = np.array(pop.tables.nodes, copy=False)
    for md in pop.diploid_metadata:
        for i in [0, 1]:
            if md.deme != nodes["deme"][md.nodes[i]]:
                return False
    return True


IndividualLocation = namedtuple("IndividualLocation", ["generation", "label", "deme"])


class WhoWhereWhen(object):
    def __init__(self):
        self.data = []

    def __call__(self, pop, sampler):
        if pop.N != len(pop.diploids):
            raise RuntimeError("incorrect total population size")
        for i in range(len(pop.diploids)):
            if i != pop.diploid_metadata[i].label:
                raise RuntimeError("index does not equal metadata label value")
            self.data.append(
                IndividualLocation(pop.generation, i, pop.diploid_metadata[i].deme)
            )


class TestConstantPopulationSize(unittest.TestCase):
    """
    This tests the absence of demographic events.
    """

    @classmethod
    def setUp(self):
        from test_demographic_models import setup_pop_rng
        from test_demographic_models import setup_pdict

        d = fwdpy11.DiscreteDemography()
        self.pop, self.rng = setup_pop_rng()
        self.pdict = setup_pdict(d, 10)
        self.params = fwdpy11.ModelParams(**self.pdict)

    def test_run_sim(self):
        fwdpy11.evolvets(self.rng, self.pop, self.params, 10)
        N = self.pop.N
        self.assertEqual(self.pop.generation, self.params.simlen)
        self.assertEqual(N, self.pop.N)


class TestSimpleMovesAndCopies(unittest.TestCase):
    """
    Basically doing the same tests
    as test_discrete_demography.TestDiscreteDemography
    """

    @classmethod
    def setUp(self):
        self.rng = fwdpy11.GSLrng(644611)
        self.pop = fwdpy11.DiploidPopulation(100, 1.0)
        self.pdict = {
            "nregions": [],
            "sregions": [],
            "recregions": [],
            "rates": (0, 0, 0),
            "demography": None,
            "simlen": 1,
            "gvalue": fwdpy11.Additive(2.0),
        }

    def test_simple_moves_from_single_deme(self):
        d = fwdpy11.DiscreteDemography([fwdpy11.move_individuals(1, 0, 1, 0.5)])
        self.pdict["demography"] = d
        self.pdict["simlen"] = 2
        params = fwdpy11.ModelParams(**self.pdict)
        w = WhoWhereWhen()
        fwdpy11.evolvets(self.rng, self.pop, params, 100, w)
        deme_sizes = self.pop.deme_sizes()
        self.assertEqual(len(deme_sizes[0]), 2)
        for i in range(len(deme_sizes[0])):
            self.assertEqual(deme_sizes[1][i], 50)
        self.assertTrue(validate_alive_node_metadata(self.pop))

        # Validate the parents (inefficiently...)
        pdata = [i for i in w.data if i.generation == self.pop.generation - 1]
        self.assertEqual(len(set(pdata)), self.pop.N)
        parents_deme_0 = set()
        parents_deme_1 = set()

        for md in self.pop.diploid_metadata:
            for p in md.parents:
                pd = [i.label for i in pdata if i.label == p]
                self.assertEqual(len(pd), 1)
                if md.deme == 0 and pd[0] not in parents_deme_0:
                    parents_deme_0.add(pd[0])
                elif md.deme == 1 and pd[0] not in parents_deme_1:
                    parents_deme_1.add(pd[0])
        self.assertEqual(len(parents_deme_0.intersection(parents_deme_1)), 0)
        self.assertEqual(len(parents_deme_1.intersection(parents_deme_0)), 0)

    def test_simple_copies_from_single_deme(self):
        d = fwdpy11.DiscreteDemography([fwdpy11.copy_individuals(1, 0, 1, 0.5)])
        self.pdict["simlen"] = 2
        self.pdict["demography"] = d
        params = fwdpy11.ModelParams(**self.pdict)
        w = WhoWhereWhen()
        fwdpy11.evolvets(self.rng, self.pop, params, 100, w)
        deme_sizes = self.pop.deme_sizes(as_dict=True)
        expected = {0: 100, 1: 50}
        self.assertEqual(deme_sizes, expected)
        self.assertTrue(validate_alive_node_metadata(self.pop))

        # Validate the parents (inefficiently...)
        pdata = [i for i in w.data if i.generation == self.pop.generation - 1]
        self.assertEqual(len(set([i.label for i in pdata])), 100)
        parents_deme_0 = set()
        parents_deme_1 = set()

        for md in self.pop.diploid_metadata:
            for p in md.parents:
                pd = [i.label for i in pdata if i.label == p]
                self.assertEqual(len(pd), 1)
                if md.deme == 0 and pd[0] not in parents_deme_0:
                    parents_deme_0.add(pd[0])
                elif md.deme == 1 and pd[0] not in parents_deme_1:
                    parents_deme_1.add(pd[0])
        self.assertGreater(len(parents_deme_0.intersection(parents_deme_1)), 0)
        self.assertGreater(len(parents_deme_1.intersection(parents_deme_0)), 0)

    def test_simple_back_and_forth_move(self):
        d = fwdpy11.DiscreteDemography(
            [
                fwdpy11.move_individuals(1, 0, 1, 0.5),
                fwdpy11.move_individuals(2, 1, 0, 1.0),
            ]
        )

        self.pdict["demography"] = d
        self.pdict["simlen"] = 3
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)
        deme_sizes = self.pop.deme_sizes()
        self.assertEqual(len(deme_sizes[0]), 1)
        for i in range(len(deme_sizes[0])):
            self.assertEqual(deme_sizes[1][i], self.pop.N)
        self.assertTrue(validate_alive_node_metadata(self.pop))

    def test_simple_back_and_forth_copy(self):
        d = fwdpy11.DiscreteDemography(
            [
                fwdpy11.copy_individuals(1, 0, 1, 0.5),
                fwdpy11.copy_individuals(2, 1, 0, 1.0),
            ]
        )

        self.pdict["demography"] = d
        self.pdict["simlen"] = 3
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)
        deme_sizes = self.pop.deme_sizes(as_dict=True)
        expected = {0: 150, 1: 50}
        self.assertEqual(deme_sizes, expected)
        self.assertTrue(validate_alive_node_metadata(self.pop))

    def test_simple_moves_from_multiple_demes(self):
        # Set 1/2 the population to start in deme 1:
        self.pop = fwdpy11.DiploidPopulation([50, 50], 1.0)
        # In generation 0, we move 1/2 of deme 1 to
        # deme 2, which creates a new deme:
        d = fwdpy11.DiscreteDemography(
            mass_migrations=[fwdpy11.move_individuals(0, 1, 2, 0.5)]
        )
        self.pdict["demography"] = d
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)
        # We now expect 50, 25, and 25 individuals in
        # demes 0, 1, and 2
        deme_sizes = self.pop.deme_sizes()
        self.assertEqual(len(deme_sizes[0]), 3)
        expected = {0: 50, 1: 25, 2: 25}
        for i in range(len(deme_sizes[0])):
            self.assertEqual(deme_sizes[1][i], expected[i])
        self.assertTrue(validate_alive_node_metadata(self.pop))

    def test_simple_copies_from_multiple_demes(self):
        # Set 1/2 the population to start in deme 1:
        self.pop = fwdpy11.DiploidPopulation([50, 50], 1.0)
        # In generation 0, we move 1/2 of deme 1 to
        # deme 2, which creates a new deme:
        d = fwdpy11.DiscreteDemography([fwdpy11.copy_individuals(0, 1, 2, 0.5)])
        self.pdict["demography"] = d
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)
        # We now expect 50, 50, and 25 individuals in
        # demes 0, 1, and 2
        deme_sizes = self.pop.deme_sizes()
        expected = {0: 50, 1: 50, 2: 25}
        self.assertEqual(len(deme_sizes[0]), len(expected))
        for i in range(len(deme_sizes[0])):
            self.assertEqual(deme_sizes[1][i], expected[i])
        self.assertTrue(validate_alive_node_metadata(self.pop))

    def test_single_deme_growth(self):
        N1 = 3412
        t = 111
        G = np.exp((np.log(N1) - np.log(self.pop.N)) / t)
        g = [fwdpy11.SetExponentialGrowth(16, 0, G)]
        d = fwdpy11.DiscreteDemography(set_growth_rates=g)
        self.pdict["demography"] = d
        self.pdict["simlen"] = 15 + t + 1
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)
        self.assertEqual(self.pop.N, N1)
        self.assertEqual(len(self.pop.diploid_metadata), N1)
        self.assertTrue(validate_alive_node_metadata(self.pop))

    def test_two_deme_growth(self):
        N0 = [90, 10]
        self.pop = fwdpy11.DiploidPopulation(N0, 1.0)
        t = [14, 23]  # generations of growth in each deme
        N1 = [5361, 616]
        G0 = np.exp((np.log(N1[0]) - np.log(N0[0])) / t[0])
        G1 = np.exp((np.log(N1[1]) - np.log(N0[1])) / t[1])

        g = []
        g.append(fwdpy11.SetExponentialGrowth(7, 0, G0))
        g.append(fwdpy11.SetExponentialGrowth(7 + t[0], 0, fwdpy11.NOGROWTH))
        g.append(fwdpy11.SetExponentialGrowth(33, 1, G1))
        g.append(fwdpy11.SetExponentialGrowth(33 + t[1], 1, fwdpy11.NOGROWTH))
        d = fwdpy11.DiscreteDemography(set_growth_rates=g)
        self.pdict["demography"] = d
        self.pdict["simlen"] = 100
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)

        md = np.array(self.pop.diploid_metadata, copy=False)
        self.assertEqual(self.pop.N, sum(N1))
        deme_sizes = self.pop.deme_sizes()
        for i, j in zip(deme_sizes[1], N1):
            self.assertEqual(i, j)
        self.assertTrue(validate_alive_node_metadata(self.pop))

    def test_two_deme_growth_with_hard_reset(self):
        N0 = [90, 10]
        self.pop = fwdpy11.DiploidPopulation(N0, 1.0)
        t = [14, 23]  # generations of growth in each deme
        N1 = [5361, 616]
        G0 = np.exp((np.log(N1[0]) - np.log(N0[0])) / t[0])
        G1 = np.exp((np.log(N1[1]) - np.log(N0[1])) / t[1])

        g = []
        g.append(fwdpy11.SetExponentialGrowth(7, 0, G0))
        g.append(fwdpy11.SetExponentialGrowth(33, 1, G1))
        g.append(fwdpy11.SetExponentialGrowth(33 + t[1], 1, fwdpy11.NOGROWTH))
        # Cut off the growth in deme 0 after a few generations,
        # and manually set the new deme size to 100 w/no growth
        p = [fwdpy11.SetDemeSize(11, 0, 100)]
        d = fwdpy11.DiscreteDemography(set_deme_sizes=p, set_growth_rates=g)
        self.pdict["demography"] = d
        self.pdict["simlen"] = 100
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)

        deme_sizes = self.pop.deme_sizes()
        N1 = [100, N1[1]]
        self.assertEqual(self.pop.N, sum(N1))
        for i, j in zip(deme_sizes[1], N1):
            self.assertEqual(i, j)
        self.assertTrue(validate_alive_node_metadata(self.pop))

    def test_two_deme_growth_without_hard_reset(self):
        N0 = [90, 10]
        self.pop = fwdpy11.DiploidPopulation(N0, 1.0)
        t = [14, 23]  # generations of growth in each deme
        N1 = [5361, 616]
        G0 = np.exp((np.log(N1[0]) - np.log(N0[0])) / t[0])
        G1 = np.exp((np.log(N1[1]) - np.log(N0[1])) / t[1])

        g = []
        g.append(fwdpy11.SetExponentialGrowth(7, 0, G0))
        g.append(fwdpy11.SetExponentialGrowth(7 + t[0], 0, fwdpy11.NOGROWTH))
        g.append(fwdpy11.SetExponentialGrowth(33, 1, G1))
        g.append(fwdpy11.SetExponentialGrowth(33 + t[1], 1, fwdpy11.NOGROWTH))
        # after X generations of growth, N[0] changes to 100
        # and the growth rate is not reset.
        p = [fwdpy11.SetDemeSize(11, 0, 100, False)]
        d = fwdpy11.DiscreteDemography(set_growth_rates=g, set_deme_sizes=p)
        self.pdict["demography"] = d
        self.pdict["simlen"] = 100
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)

        md = np.array(self.pop.diploid_metadata, copy=False)
        deme_sizes = self.pop.deme_sizes()
        N1[0] = np.round(100.0 * np.power(G0, 7 + t[0] - 11))
        self.assertEqual(self.pop.N, sum(N1))
        for i, j in zip(deme_sizes[1], N1):
            self.assertEqual(i, j)
        self.assertTrue(validate_alive_node_metadata(self.pop))

    def test_two_moves_in_same_generation(self):
        """
        In generation 0, deme 0 splits equally
        into demes 1 and 2.  This should leave
        deme 0 empty and the sizes of demes
        1 and 2 both equal to 0.5 the initial
        size.
        """
        m = [
            fwdpy11.move_individuals(0, 0, 2, 0.5),
            fwdpy11.move_individuals(0, 0, 1, 0.5),
        ]
        d = fwdpy11.DiscreteDemography(m)
        N0 = self.pop.N
        self.pdict["demography"] = d
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)
        expected = {0: 0, 1: N0 // 2, 2: N0 // 2}
        deme_sizes = self.pop.deme_sizes()
        for i, j in zip(deme_sizes[0], deme_sizes[1]):
            self.assertEqual(expected[i], j)
        self.assertTrue(validate_alive_node_metadata(self.pop))

    def test_copies_happen_before_moves(self):
        """
        This is a more complex test.
        In the same generation:
        1. All individuals are copied from deme 0 to deme 1.
        2. 50% of deme 0 moves to deme 2.
        3. 50% of deme 0 moves to deme 3.

        This, at the end, the size of deme 1 should be
        the initial size and the size of demes 2 and 3
        should be 0.5*initial size
        """
        m = [
            fwdpy11.move_individuals(0, 0, 3, 0.5),
            fwdpy11.copy_individuals(0, 0, 1, 1),
            fwdpy11.move_individuals(0, 0, 2, 0.5),
        ]
        d = fwdpy11.DiscreteDemography(m)
        N0 = self.pop.N
        self.pdict["demography"] = d
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)
        expected = {0: 0, 1: N0, 2: N0 // 2, 3: N0 // 2}
        deme_sizes = self.pop.deme_sizes()
        for i, j in zip(deme_sizes[0], deme_sizes[1]):
            self.assertEqual(expected[i], j)
        self.assertTrue(validate_alive_node_metadata(self.pop))

    def test_mass_move_with_growth(self):
        """
        In generation 5, the mass movement from 0 to 1
        will reset the growth rate in deme 0 to
        fwdpy11.NOGROWTH

        This test is also handy b/c growth results in
        an odd total N by the time the mass migration happens
        """
        m = [fwdpy11.move_individuals(5, 0, 1, 0.5)]
        g = [fwdpy11.SetExponentialGrowth(0, 0, 1.1)]
        d = fwdpy11.DiscreteDemography(mass_migrations=m, set_growth_rates=g)
        N0 = self.pop.N
        self.pdict["demography"] = d
        self.pdict["simlen"] = 10
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)
        deme_sizes = self.pop.deme_sizes()
        N5 = np.round(N0 * np.power(g[0].G, 5))
        N_after_mass_mig = [N5 // 2, N5 - N5 // 2]
        for i, j in zip(deme_sizes[1], N_after_mass_mig):
            self.assertEqual(i, j)
        self.assertTrue(validate_alive_node_metadata(self.pop))

    def test_mass_move_with_growth_no_reset(self):
        """
        In generation 5, the mass movement from 0 to 1
        """
        m = [fwdpy11.move_individuals(5, 0, 1, 0.5, False)]
        g = [fwdpy11.SetExponentialGrowth(0, 0, 1.1)]
        d = fwdpy11.DiscreteDemography(mass_migrations=m, set_growth_rates=g)
        N0 = self.pop.N
        self.pdict["demography"] = d
        self.pdict["simlen"] = 10
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)
        deme_sizes = self.pop.deme_sizes()
        N5 = np.round(N0 * np.power(g[0].G, 5))
        N_after_mass_mig_0 = np.round((N5 // 2) * np.power(g[0].G, 5))
        N = [N_after_mass_mig_0, N5 - N5 // 2]
        for i, j in zip(N, deme_sizes[1]):
            self.assertEqual(i, j)
        self.assertTrue(validate_alive_node_metadata(self.pop))


class TestSimpleDemeSizeChanges(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.rng = fwdpy11.GSLrng(42)
        self.pop = fwdpy11.DiploidPopulation(100, 1.0)
        self.pdict = {
            "nregions": [],
            "sregions": [],
            "recregions": [],
            "rates": (0, 0, 0),
            "demography": None,
            "simlen": 1,
            "gvalue": fwdpy11.Additive(2.0),
        }

    def test_global_extinction_single_deme(self):
        ds = [fwdpy11.SetDemeSize(4, 0, 0)]
        d = fwdpy11.DiscreteDemography(set_deme_sizes=ds)
        self.pdict["demography"] = d
        self.pdict["simlen"] = 10
        params = fwdpy11.ModelParams(**self.pdict)
        with self.assertRaises(fwdpy11.GlobalExtinction):
            fwdpy11.evolvets(self.rng, self.pop, params, 100)

        self.assertEqual(self.pop.N, 100)
        self.assertEqual(self.pop.generation, 4)
        nodes = np.array(self.pop.tables.nodes)
        self.assertTrue(np.all(nodes["time"][self.pop.alive_nodes] == 4))

    def test_no_valid_parents(self):
        m = [fwdpy11.move_individuals(0, 0, 1, 1)]
        s = [fwdpy11.SetDemeSize(0, 0, 100)]
        d = fwdpy11.DiscreteDemography(mass_migrations=m, set_deme_sizes=s)
        self.pdict["demography"] = d
        self.pdict["simlen"] = 10
        params = fwdpy11.ModelParams(**self.pdict)
        with self.assertRaises(fwdpy11.DemographyError):
            fwdpy11.evolvets(self.rng, self.pop, params, 100)

    def test_no_valid_parents_alternate_method(self):
        m = [fwdpy11.copy_individuals(0, 0, 1, 1.0)]
        s = [fwdpy11.SetDemeSize(0, 0, 0), fwdpy11.SetDemeSize(2, 0, 100)]
        d = fwdpy11.DiscreteDemography(mass_migrations=m, set_deme_sizes=s)
        self.pdict["demography"] = d
        self.pdict["simlen"] = 10
        params = fwdpy11.ModelParams(**self.pdict)
        with self.assertRaises(fwdpy11.DemographyError):
            fwdpy11.evolvets(self.rng, self.pop, params, 100)

    def test_no_valid_parents_with_migration(self):
        m = [fwdpy11.move_individuals(0, 0, 1, 1)]
        s = [fwdpy11.SetDemeSize(0, 0, 100)]
        M = np.array([0.5] * 4).reshape(2, 2)
        d = fwdpy11.DiscreteDemography(mass_migrations=m, set_deme_sizes=s, migmatrix=M)
        self.pdict["demography"] = d
        self.pdict["simlen"] = 10
        params = fwdpy11.ModelParams(**self.pdict)
        with self.assertRaises(fwdpy11.DemographyError):
            fwdpy11.evolvets(self.rng, self.pop, params, 100)

    def test_no_valid_parents_with_migration_2(self):
        m = [fwdpy11.move_individuals(0, 0, 1, 1)]
        s = [fwdpy11.SetDemeSize(0, 0, 100)]
        M = np.array([0.5] * 4).reshape(2, 2)
        cM = [fwdpy11.SetMigrationRates(0, None, np.array([0, 1, 0, 1]).reshape(2, 2))]
        d = fwdpy11.DiscreteDemography(
            mass_migrations=m, set_deme_sizes=s, set_migration_rates=cM, migmatrix=M
        )
        self.pdict["demography"] = d
        self.pdict["simlen"] = 10
        params = fwdpy11.ModelParams(**self.pdict)
        try:
            fwdpy11.evolvets(self.rng, self.pop, params, 100)
        except Exception:
            self.fail("unexpected exception")


class TestSimpleMigrationModels(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.rng = fwdpy11.GSLrng(42)
        self.pop = fwdpy11.DiploidPopulation(100, 1.0)
        self.pdict = {
            "nregions": [],
            "sregions": [],
            "recregions": [],
            "rates": (0, 0, 0),
            "demography": None,
            "simlen": 1,
            "gvalue": fwdpy11.Additive(2.0),
        }

    def test_simple_two_deme_migration(self):
        mm = np.array([0.5] * 4).reshape(2, 2)
        self.pop = fwdpy11.DiploidPopulation([50, 50], 1.0)
        d = fwdpy11.DiscreteDemography(migmatrix=mm)
        N = self.pop.N
        self.pdict["demography"] = d
        self.pdict["simlen"] = 5
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)
        self.assertEqual(N, self.pop.N)
        deme_sizes = self.pop.deme_sizes()
        for i in deme_sizes[1]:
            self.assertEqual(i, self.pop.N // 2)
        self.assertTrue(validate_alive_node_metadata(self.pop))

    def test_change_migration_rates_simple_two_deme_migration(self):
        """
        For a 2-deme model, the mig matrix is
        [[0, 1]
         [1, 0]]
        so that all offspring have both parents from the other deme.

        After 3 generations, we reset the migration rates to be
        [[1, 0],
         [1, 0]]
        so that all parents are from deme zero.
        """
        mm = np.array([0, 1, 1, 0]).reshape(2, 2)
        mmigs = [fwdpy11.move_individuals(0, 0, 1, 0.5)]
        smr = [fwdpy11.SetMigrationRates(3, None, np.array([1, 0, 1, 0]).reshape(2, 2))]
        d = fwdpy11.DiscreteDemography(
            mass_migrations=mmigs, migmatrix=mm, set_migration_rates=smr
        )
        N = self.pop.N
        self.pdict["demography"] = d
        self.pdict["simlen"] = 5
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)
        self.assertEqual(N, self.pop.N)
        deme_sizes = self.pop.deme_sizes()
        for i in deme_sizes[1]:
            self.assertEqual(i, self.pop.N // 2)
        self.assertTrue(validate_alive_node_metadata(self.pop))

    def test_change_migration_rates_simple_two_deme_migration_bad_matrix(self):
        """
        For a 2-deme model, the mig matrix is
        [0, 1
         1, 0]
        so that all offspring have both parents from the other deme.

        After 3 generations, we reset the migration rates to be
        [[0.5, 0.5],
         [0, 0]],
        which leads to there being no parents for deme 1, raising a
        fwdpy11.DemographyError exception.
        """
        mm = np.array([0, 1, 1, 0]).reshape(2, 2)
        mmigs = [fwdpy11.move_individuals(0, 0, 1, 0.5)]
        smr = [
            fwdpy11.SetMigrationRates(3, None, np.array([0.5, 0.5, 0, 0]).reshape(2, 2))
        ]
        d = fwdpy11.DiscreteDemography(
            mass_migrations=mmigs, migmatrix=mm, set_migration_rates=smr
        )
        self.pdict["demography"] = d
        self.pdict["simlen"] = 5
        params = fwdpy11.ModelParams(**self.pdict)
        with self.assertRaises(fwdpy11.DemographyError):
            fwdpy11.evolvets(self.rng, self.pop, params, 100)

    def test_selfing_vs_migration(self):
        """
        Parents of deme 0 are all migrants from deme 1.
        Parents of deme 1 are all migrants from deme 0.
        Deme 0 never selfs.  Deme 1 always selfs.

        NOTE: currently doesn't test anything.
        """
        mm = np.array([0, 1, 1, 0]).reshape(2, 2)
        mmigs = [fwdpy11.move_individuals(0, 0, 1, 0.5)]
        s = [fwdpy11.SetSelfingRate(0, 1, 1.0)]
        mmigs = [fwdpy11.move_individuals(0, 0, 1, 0.5)]
        d = fwdpy11.DiscreteDemography(
            mass_migrations=mmigs, migmatrix=mm, set_selfing_rates=s
        )
        self.pdict["demography"] = d
        self.pdict["simlen"] = 5
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)
        self.assertTrue(validate_alive_node_metadata(self.pop))


class TestGrowthModels(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.rng = fwdpy11.GSLrng(42)
        self.pop = fwdpy11.DiploidPopulation(100, 1.0)
        self.pdict = {
            "nregions": [],
            "sregions": [],
            "recregions": [],
            "rates": (0, 0, 0),
            "demography": None,
            "simlen": 1,
            "gvalue": fwdpy11.Additive(2.0),
        }

    def test_global_extinction_single_deme(self):
        """
        Pop would go exctinct in generation 8,
        which triggers exception in generation 7
        """
        g = [fwdpy11.SetExponentialGrowth(0, 0, 0.5)]
        d = fwdpy11.DiscreteDemography(set_growth_rates=g)
        self.pdict["demography"] = d
        self.pdict["simlen"] = 20
        params = fwdpy11.ModelParams(**self.pdict)
        with self.assertRaises(fwdpy11.GlobalExtinction):
            fwdpy11.evolvets(self.rng, self.pop, params, 100)
        self.assertEqual(self.pop.generation, 7)

        nodes = np.array(self.pop.tables.nodes)
        self.assertTrue(np.all(nodes["time"][self.pop.alive_nodes] == 7))

    def test_shrink_to_zero_then_recolonize(self):
        m = [
            fwdpy11.copy_individuals(0, 0, 1, 1.0),
            fwdpy11.copy_individuals(10, 1, 0, 1.0),
        ]
        g = [fwdpy11.SetExponentialGrowth(0, 0, 0.5)]
        d = fwdpy11.DiscreteDemography(mass_migrations=m, set_growth_rates=g)
        self.pdict["demography"] = d
        self.pdict["simlen"] = 20
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)
        deme_sizes = self.pop.deme_sizes(as_dict=True)
        self.assertEqual(deme_sizes, {0: 100, 1: 100})


class TestGeneticValueLists(unittest.TestCase):
    """
    Added in 0.6.0 to test exceptions raised
    when the number of genetic value objects
    are out of whack vis-a-vis the demography.
    """

    @classmethod
    def setUp(self):
        self.rng = fwdpy11.GSLrng(42)
        self.pop = fwdpy11.DiploidPopulation(100, 1.0)
        self.pdict = {
            "nregions": [],
            "sregions": [],
            "recregions": [],
            "rates": (0, 0, 0),
            "demography": fwdpy11.DiscreteDemography(),
            "simlen": 1,
            "gvalue": fwdpy11.Additive(2.0),
        }

    def test_too_many(self):
        self.pdict["gvalue"] = [fwdpy11.Additive(2.0)] * 2
        params = fwdpy11.ModelParams(**self.pdict)
        with self.assertRaises(ValueError):
            fwdpy11.evolvets(self.rng, self.pop, params, 100)

    def test_too_few(self):
        mmigs = [
            fwdpy11.move_individuals(0, 0, 1, 1.0 / 3.0),
            fwdpy11.move_individuals(0, 0, 2, 1.0 / 3.0),
        ]
        d = fwdpy11.DiscreteDemography(mass_migrations=mmigs)
        self.pdict["demography"] = d
        self.pdict["simlen"] = 5
        self.pdict["gvalue"] = [fwdpy11.Additive(2.0)] * 2  # This is the error
        params = fwdpy11.ModelParams(**self.pdict)
        with self.assertRaises(ValueError):
            fwdpy11.evolvets(self.rng, self.pop, params, 100)

    def test_too_few_alt_method(self):
        self.pop = fwdpy11.DiploidPopulation(
            [0, self.pop.N // 3, 2 * self.pop.N // 3], 1
        )
        self.pdict["gvalue"] = [fwdpy11.Additive(2.0)] * 2  # This is the error
        params = fwdpy11.ModelParams(**self.pdict)
        with self.assertRaises(ValueError):
            fwdpy11.evolvets(self.rng, self.pop, params, 100)


class TestIMModel(unittest.TestCase):
    """
    Test of an Isolation-with-Migration, or IM model.

    We also test that we can evolve a pop up to a certain
    point, stop, then evolve it "the rest of the way". Such
    a test is useful for simulations of human OOA models
    where we may want to simplify more often during the
    exponential growth phase in order to save memory.
    """

    def getG(N0, Nt, t):
        return np.exp((np.log(Nt) - np.log(N0)) / t)

    @classmethod
    def setUp(self):
        self.rng = fwdpy11.GSLrng(666)
        self.Nref = 200
        self.Tsplit = 10
        self.psplit = 0.2
        self.gens_post_split = 5
        self.N0t = 300
        self.N1t = 400
        G0 = self.getG(self.Nref * (1.0 - self.psplit), self.N0t, self.gens_post_split)
        G1 = self.getG(self.Nref * (self.psplit), self.N1t, self.gens_post_split)

        split = [fwdpy11.move_individuals(self.Tsplit, 0, 1, self.psplit)]
        growth = [
            fwdpy11.SetExponentialGrowth(self.Tsplit, 0, G0),
            fwdpy11.SetExponentialGrowth(self.Tsplit, 1, G1),
        ]
        m = np.zeros(4).reshape(2, 2)
        m[0, 0] = 1
        cm = [
            fwdpy11.SetMigrationRates(self.Tsplit, 0, [0.98, 0.02]),
            fwdpy11.SetMigrationRates(self.Tsplit, 1, [0.02, 0.98]),
        ]
        d = fwdpy11.DiscreteDemography(
            mass_migrations=split,
            set_growth_rates=growth,
            set_migration_rates=cm,
            migmatrix=m,
        )

        nregions = []
        sregions = []
        recregions = []

        self.pdict = {
            "nregions": nregions,
            "sregions": sregions,
            "recregions": recregions,
            "rates": (0, 0, 0),
            "gvalue": fwdpy11.Multiplicative(2.0),
            "demography": d,
            "simlen": -1,
            "prune_selected": True,
        }
        self.pop = fwdpy11.DiploidPopulation(self.Nref, 1.0)

    def test_complete_sim(self):
        self.pdict["simlen"] = self.Tsplit + self.gens_post_split
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 5)
        self.assertEqual(self.pop.generation, params.simlen)
        deme_sizes = self.pop.deme_sizes()
        self.assertEqual(deme_sizes[1][0], self.N0t)
        self.assertEqual(deme_sizes[1][1], self.N1t)

    def test_evolve_in_two_steps(self):
        self.pdict["simlen"] = self.Tsplit
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 5)
        self.assertEqual(self.pop.generation, self.Tsplit)
        deme_sizes = self.pop.deme_sizes()
        self.assertEqual(len(deme_sizes[0]), 2)
        self.assertEqual(deme_sizes[1][0], int((1.0 - self.psplit) * self.Nref))

        self.pdict["simlen"] = self.gens_post_split
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 2)  # Simplify more often
        self.assertEqual(self.pop.generation, self.Tsplit + self.gens_post_split)
        deme_sizes = self.pop.deme_sizes()
        self.assertEqual(deme_sizes[1][0], self.N0t)
        self.assertEqual(deme_sizes[1][1], self.N1t)

    def test_evolve_in_two_steps_restart_with_two_demes(self):
        """
        Tests the more complex case of restarting a sim
        when multiple demes are present.
        """
        self.pdict["simlen"] = self.Tsplit + 3  # Evolve past the split
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 5)
        self.assertEqual(self.pop.generation, self.Tsplit + 3)
        deme_sizes = self.pop.deme_sizes()
        self.assertTrue(all([i > 0 for i in deme_sizes[1].tolist()]))

        self.pdict["simlen"] = self.gens_post_split - 3
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(
            self.rng, self.pop, params, 2, check_demographic_event_timings=False
        )  # Simplify more often
        self.assertEqual(self.pop.generation, self.Tsplit + self.gens_post_split)
        deme_sizes = self.pop.deme_sizes()
        self.assertEqual(deme_sizes[1][0], self.N0t)
        self.assertEqual(deme_sizes[1][1], self.N1t)

    def test_evolve_in_two_steps_restart_with_two_demes_and_pickle(self):
        """
        Tests the more complex case of restarting a sim
        when multiple demes are present and we pickle the
        demographic model.
        """
        self.pdict["simlen"] = self.Tsplit + 3  # Evolve past the split
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 5)
        self.assertEqual(self.pop.generation, self.Tsplit + 3)
        deme_sizes = self.pop.deme_sizes()
        self.assertTrue(all([i > 0 for i in deme_sizes[1].tolist()]))

        self.pdict["simlen"] = self.gens_post_split - 3
        params = fwdpy11.ModelParams(**self.pdict)
        pparams = pickle.dumps(params, -1)
        unpickled_params = pickle.loads(pparams)
        fwdpy11.evolvets(
            self.rng,
            self.pop,
            unpickled_params,
            2,
            check_demographic_event_timings=False,
        )  # Simplify more often
        self.assertEqual(self.pop.generation, self.Tsplit + self.gens_post_split)
        deme_sizes = self.pop.deme_sizes()
        self.assertEqual(deme_sizes[1][0], self.N0t)
        self.assertEqual(deme_sizes[1][1], self.N1t)

    def test_reuse_of_the_demographic_model(self):
        """
        Since DiscreteDemography "remembers" the state
        of the model, we have to test that the model gets
        properly reset when confronted with a "fresh"
        population.
        """
        self.pdict["simlen"] = self.Tsplit + self.gens_post_split
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 5)
        self.assertEqual(self.pop.generation, params.simlen)
        deme_sizes = self.pop.deme_sizes()
        self.assertEqual(deme_sizes[1][0], self.N0t)
        self.assertEqual(deme_sizes[1][1], self.N1t)

        npop = fwdpy11.DiploidPopulation(self.Nref, 1.0)
        fwdpy11.evolvets(self.rng, npop, params, 5)
        self.assertEqual(npop.generation, params.simlen)
        deme_sizes = self.pop.deme_sizes()
        self.assertEqual(deme_sizes[1][0], self.N0t)
        self.assertEqual(deme_sizes[1][1], self.N1t)


if __name__ == "__main__":
    unittest.main()
