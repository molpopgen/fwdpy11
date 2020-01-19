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
import unittest
import fwdpy11
import numpy as np


def validate_alive_node_metadata(pop):
    nodes = np.array(pop.tables.nodes, copy=False)
    for md in pop.diploid_metadata:
        for i in [0, 1]:
            if md.deme != nodes['deme'][md.nodes[i]]:
                return False
    return True


class TestSimpleMovesAndCopies(unittest.TestCase):
    """
    Basically doing the same tests
    as test_discrete_demography.TestDiscreteDemography
    """
    @classmethod
    def setUp(self):
        self.rng = fwdpy11.GSLrng(42)
        self.pop = fwdpy11.DiploidPopulation(100, 1.)
        self.pdict = {'nregions': [],
                      'sregions': [],
                      'recregions': [],
                      'rates': (0, 0, 0),
                      'demography': None,
                      'simlen': 1,
                      'gvalue': fwdpy11.Additive(2.)
                      }

    def test_simple_moves_from_single_deme(self):
        d = fwdpy11.DiscreteDemography(
            [fwdpy11.move_individuals(0, 0, 1, 0.5)])
        self.pdict['demography'] = d
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)
        md = np.array(self.pop.diploid_metadata, copy=False)
        deme_counts = np.unique(md['deme'], return_counts=True)
        self.assertEqual(len(deme_counts[0]), 2)
        for i in range(len(deme_counts[0])):
            self.assertEqual(deme_counts[1][i], 50)
        self.assertTrue(validate_alive_node_metadata(self.pop))

    def test_simple_copies_from_single_deme(self):
        d = fwdpy11.DiscreteDemography(
            [fwdpy11.copy_individuals(0, 0, 1, 0.5)])
        self.pdict['demography'] = d
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)
        md = np.array(self.pop.diploid_metadata, copy=False)
        deme_counts = np.unique(md['deme'], return_counts=True)
        expected = {0: 100, 1: 50}
        self.assertEqual(len(deme_counts[0]), len(expected))
        for i in range(len(deme_counts[0])):
            self.assertEqual(deme_counts[1][i], expected[i])
        self.assertTrue(validate_alive_node_metadata(self.pop))

    def test_simple_back_and_forth_move(self):
        d = fwdpy11.DiscreteDemography(
            [fwdpy11.move_individuals(1, 0, 1, 0.5),
             fwdpy11.move_individuals(2, 1, 0, 1.)])

        self.pdict['demography'] = d
        self.pdict['simlen'] = 3
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)
        md = np.array(self.pop.diploid_metadata, copy=False)
        deme_counts = np.unique(md['deme'], return_counts=True)
        self.assertEqual(len(deme_counts[0]), 1)
        for i in range(len(deme_counts[0])):
            self.assertEqual(deme_counts[1][i], self.pop.N)
        self.assertTrue(validate_alive_node_metadata(self.pop))

    def test_simple_back_and_forth_copy(self):
        d = fwdpy11.DiscreteDemography(
            [fwdpy11.copy_individuals(1, 0, 1, 0.5),
             fwdpy11.copy_individuals(2, 1, 0, 1.)])

        self.pdict['demography'] = d
        self.pdict['simlen'] = 3
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)
        md = np.array(self.pop.diploid_metadata, copy=False)
        deme_counts = np.unique(md['deme'], return_counts=True)
        expected = {0: 150, 1: 50}
        self.assertEqual(len(deme_counts[0]), len(expected))
        for i in range(len(deme_counts[0])):
            self.assertEqual(deme_counts[1][i], expected[i])
        self.assertTrue(validate_alive_node_metadata(self.pop))

    def test_simple_moves_from_multiple_demes(self):
        # Set 1/2 the population to start in deme 1:
        md = np.array(self.pop.diploid_metadata, copy=False)
        md['deme'][self.pop.N//2:] = 1
        # In generation 0, we move 1/2 of deme 1 to
        # deme 2, which creates a new deme:
        d = fwdpy11.DiscreteDemography(
            [fwdpy11.move_individuals(0, 1, 2, 0.5)])
        self.pdict['demography'] = d
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)
        # We now expect 50, 25, and 25 individuals in
        # demes 0, 1, and 2
        md = np.array(self.pop.diploid_metadata, copy=False)
        deme_counts = np.unique(md['deme'], return_counts=True)
        self.assertEqual(len(deme_counts[0]), 3)
        expected = {0: 50, 1: 25, 2: 25}
        for i in range(len(deme_counts[0])):
            self.assertEqual(deme_counts[1][i], expected[i])
        self.assertTrue(validate_alive_node_metadata(self.pop))

    def test_simple_copies_from_multiple_demes(self):
        # Set 1/2 the population to start in deme 1:
        md = np.array(self.pop.diploid_metadata, copy=False)
        md['deme'][self.pop.N//2:] = 1
        # In generation 0, we move 1/2 of deme 1 to
        # deme 2, which creates a new deme:
        d = fwdpy11.DiscreteDemography(
            [fwdpy11.copy_individuals(0, 1, 2, 0.5)])
        self.pdict['demography'] = d
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)
        # We now expect 50, 50, and 25 individuals in
        # demes 0, 1, and 2
        md = np.array(self.pop.diploid_metadata, copy=False)
        deme_counts = np.unique(md['deme'], return_counts=True)
        expected = {0: 50, 1: 50, 2: 25}
        self.assertEqual(len(deme_counts[0]), len(expected))
        for i in range(len(deme_counts[0])):
            self.assertEqual(deme_counts[1][i], expected[i])
        self.assertTrue(validate_alive_node_metadata(self.pop))

    def test_single_deme_growth(self):
        N1 = 3412
        t = 111
        G = np.exp((np.log(N1)-np.log(self.pop.N))/t)
        g = [fwdpy11.SetExponentialGrowth(16, 0, G)]
        d = fwdpy11.DiscreteDemography(set_growth_rates=g)
        self.pdict['demography'] = d
        self.pdict['simlen'] = 15+t+1
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)
        self.assertEqual(self.pop.N, N1)
        self.assertEqual(len(self.pop.diploid_metadata), N1)
        self.assertTrue(validate_alive_node_metadata(self.pop))

    def test_two_deme_growth(self):
        N0 = [90, 10]
        md = np.array(self.pop.diploid_metadata, copy=False)
        md['deme'][N0[0]:] = 1
        t = [14, 23]  # generations of growth in each deme
        N1 = [5361, 616]
        G0 = np.exp((np.log(N1[0])-np.log(N0[0]))/t[0])
        G1 = np.exp((np.log(N1[1])-np.log(N0[1]))/t[1])

        g = []
        g.append(fwdpy11.SetExponentialGrowth(7, 0, G0))
        g.append(fwdpy11.SetExponentialGrowth(7+t[0], 0, fwdpy11.NOGROWTH))
        g.append(fwdpy11.SetExponentialGrowth(33, 1, G1))
        g.append(fwdpy11.SetExponentialGrowth(33+t[1], 1, fwdpy11.NOGROWTH))
        d = fwdpy11.DiscreteDemography(set_growth_rates=g)
        self.pdict['demography'] = d
        self.pdict['simlen'] = 100
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)

        md = np.array(self.pop.diploid_metadata, copy=False)
        self.assertEqual(self.pop.N, sum(N1))
        deme_sizes = np.unique(md['deme'], return_counts=True)
        for i, j in zip(deme_sizes[1], N1):
            self.assertEqual(i, j)
        self.assertTrue(validate_alive_node_metadata(self.pop))

    def test_two_deme_growth_with_hard_reset(self):
        N0 = [90, 10]
        md = np.array(self.pop.diploid_metadata, copy=False)
        md['deme'][N0[0]:] = 1
        t = [14, 23]  # generations of growth in each deme
        N1 = [5361, 616]
        G0 = np.exp((np.log(N1[0])-np.log(N0[0]))/t[0])
        G1 = np.exp((np.log(N1[1])-np.log(N0[1]))/t[1])

        g = []
        g.append(fwdpy11.SetExponentialGrowth(7, 0, G0))
        g.append(fwdpy11.SetExponentialGrowth(33, 1, G1))
        g.append(fwdpy11.SetExponentialGrowth(33+t[1], 1, fwdpy11.NOGROWTH))
        # Cut off the growth in deme 0 after a few generations,
        # and manually set the new deme size to 100 w/no growth
        p = [fwdpy11.SetDemeSize(11, 0, 100)]
        d = fwdpy11.DiscreteDemography(set_deme_sizes=p,
                                       set_growth_rates=g)
        self.pdict['demography'] = d
        self.pdict['simlen'] = 100
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)

        md = np.array(self.pop.diploid_metadata, copy=False)
        deme_sizes = np.unique(md['deme'], return_counts=True)
        N1 = [100, N1[1]]
        self.assertEqual(self.pop.N, sum(N1))
        for i, j in zip(deme_sizes[1], N1):
            self.assertEqual(i, j)
        self.assertTrue(validate_alive_node_metadata(self.pop))

    def test_two_deme_growth_without_hard_reset(self):
        N0 = [90, 10]
        md = np.array(self.pop.diploid_metadata, copy=False)
        md['deme'][N0[0]:] = 1
        t = [14, 23]  # generations of growth in each deme
        N1 = [5361, 616]
        G0 = np.exp((np.log(N1[0])-np.log(N0[0]))/t[0])
        G1 = np.exp((np.log(N1[1])-np.log(N0[1]))/t[1])

        g = []
        g.append(fwdpy11.SetExponentialGrowth(7, 0, G0))
        g.append(fwdpy11.SetExponentialGrowth(7+t[0], 0, fwdpy11.NOGROWTH))
        g.append(fwdpy11.SetExponentialGrowth(33, 1, G1))
        g.append(fwdpy11.SetExponentialGrowth(33+t[1], 1, fwdpy11.NOGROWTH))
        # after X generations of growth, N[0] changes to 100
        # and the growth rate is not reset.
        p = [fwdpy11.SetDemeSize(11, 0, 100, False)]
        d = fwdpy11.DiscreteDemography(set_growth_rates=g,
                                       set_deme_sizes=p)
        self.pdict['demography'] = d
        self.pdict['simlen'] = 100
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)

        md = np.array(self.pop.diploid_metadata, copy=False)
        deme_sizes = np.unique(md['deme'], return_counts=True)
        N1[0] = np.round(100.*np.power(G0, 7 + t[0] - 11))
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
        m = [fwdpy11.move_individuals(0, 0, 2, 0.5),
             fwdpy11.move_individuals(0, 0, 1, 0.5)]
        d = fwdpy11.DiscreteDemography(m)
        N0 = self.pop.N
        self.pdict['demography'] = d
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)
        md = np.array(self.pop.diploid_metadata, copy=False)
        expected = {0: 0, 1: N0//2, 2: N0//2}
        dc = np.unique(md['deme'], return_counts=True)
        for i, j in zip(dc[0], dc[1]):
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
        m = [fwdpy11.move_individuals(0, 0, 3, 0.5),
             fwdpy11.copy_individuals(0, 0, 1, 1),
             fwdpy11.move_individuals(0, 0, 2, 0.5)]
        d = fwdpy11.DiscreteDemography(m)
        N0 = self.pop.N
        self.pdict['demography'] = d
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)
        md = np.array(self.pop.diploid_metadata, copy=False)
        expected = {0: 0, 1: N0, 2: N0//2, 3: N0//2}
        dc = np.unique(md['deme'], return_counts=True)
        for i, j in zip(dc[0], dc[1]):
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
        self.pdict['demography'] = d
        self.pdict['simlen'] = 10
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)
        md = np.array(self.pop.diploid_metadata)
        deme_counts = np.unique(md['deme'], return_counts=True)
        N5 = np.round(N0*np.power(g[0].G, 5))
        N_after_mass_mig = [N5//2, N5-N5//2]
        for i, j in zip(deme_counts[1], N_after_mass_mig):
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
        self.pdict['demography'] = d
        self.pdict['simlen'] = 10
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)
        md = np.array(self.pop.diploid_metadata)
        deme_counts = np.unique(md['deme'], return_counts=True)
        N5 = np.round(N0*np.power(g[0].G, 5))
        N_after_mass_mig_0 = np.round((N5//2)*np.power(g[0].G, 5))
        N = [N_after_mass_mig_0, N5 - N5//2]
        for i, j in zip(N, deme_counts[1]):
            self.assertEqual(i, j)
        self.assertTrue(validate_alive_node_metadata(self.pop))


class TestSimpleMigrationModels(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.rng = fwdpy11.GSLrng(42)
        self.pop = fwdpy11.DiploidPopulation(100, 1.)
        self.pdict = {'nregions': [],
                      'sregions': [],
                      'recregions': [],
                      'rates': (0, 0, 0),
                      'demography': None,
                      'simlen': 1,
                      'gvalue': fwdpy11.Additive(2.)
                      }

    def test_simple_two_deme_migration(self):
        mm = np.array([0.5]*4).reshape(2, 2)
        md = np.array(self.pop.diploid_metadata, copy=False)
        md['deme'][len(md)//2:] = 1
        d = fwdpy11.DiscreteDemography(migmatrix=mm)
        N = self.pop.N
        self.pdict['demography'] = d
        self.pdict['simlen'] = 5
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)
        self.assertEqual(N, self.pop.N)
        md = np.array(self.pop.diploid_metadata, copy=False)
        deme_sizes = np.unique(md['deme'], return_counts=True)
        for i in deme_sizes[1]:
            self.assertEqual(i, self.pop.N//2)
        self.assertTrue(validate_alive_node_metadata(self.pop))

    def test_change_migration_rates_simple_two_deme_migration(self):
        """
        For a 2-deme model, the mig matrix is
        [[0, 1]
         [1, 0]]
        so that all offspring have both parents from the other deme.

        After 3 generations, we reset the migration rates to be
        [[0.5, 0.5],
         [[0, 0]]
        so that all parents are from deme zero.
        """
        mm = np.array([0, 1, 1, 0]).reshape(2, 2)
        mmigs = [fwdpy11.move_individuals(0, 0, 1, 0.5)]
        smr = [fwdpy11.SetMigrationRates(
            3, np.array([0.5, 0.5, 0, 0]).reshape(2, 2))]
        d = fwdpy11.DiscreteDemography(mass_migrations=mmigs, migmatrix=mm,
                                       set_migration_rates=smr)
        N = self.pop.N
        self.pdict['demography'] = d
        self.pdict['simlen'] = 5
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)
        self.assertEqual(N, self.pop.N)
        md = np.array(self.pop.diploid_metadata, copy=False)
        deme_sizes = np.unique(md['deme'], return_counts=True)
        for i in deme_sizes[1]:
            self.assertEqual(i, self.pop.N//2)
        self.assertTrue(validate_alive_node_metadata(self.pop))

    def test_change_migration_rates_simple_two_deme_migration_bad_matrix(self):
        """
        For a 2-deme model, the mig matrix is
        [0, 1
         1, 0]
        so that all offspring have both parents from the other deme.

        After 3 generations, we reset the migration rates to be
        [[1,0],
         [1,0]],
        which leads to there being no parents for deme 1, raising a
        MigrationError exception.
        """
        mm = np.array([0, 1, 1, 0]).reshape(2, 2)
        mmigs = [fwdpy11.move_individuals(0, 0, 1, 0.5)]
        smr = [fwdpy11.SetMigrationRates(
            3, np.array([1, 0, 1, 0]).reshape(2, 2))]
        d = fwdpy11.DiscreteDemography(mass_migrations=mmigs, migmatrix=mm,
                                       set_migration_rates=smr)
        self.pdict['demography'] = d
        self.pdict['simlen'] = 5
        params = fwdpy11.ModelParams(**self.pdict)
        with self.assertRaises(fwdpy11.MigrationError):
            fwdpy11.evolvets(self.rng, self.pop, params, 100)
        self.assertTrue(validate_alive_node_metadata(self.pop))

    def test_migrration_rates_larger_than_one(self):
        """
        Same as a previous tests, but rates are "weights"
        rather than "probabilities"

        """
        mm = np.array([0, 2, 2, 0]).reshape(2, 2)
        mmigs = [fwdpy11.move_individuals(0, 0, 1, 0.5)]
        smr = [fwdpy11.SetMigrationRates(
            3, np.array([1.5, 1.5, 0, 0]).reshape(2, 2))]
        d = fwdpy11.DiscreteDemography(mass_migrations=mmigs, migmatrix=mm,
                                       set_migration_rates=smr)
        N = self.pop.N
        self.pdict['demography'] = d
        self.pdict['simlen'] = 5
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)
        self.assertEqual(N, self.pop.N)
        md = np.array(self.pop.diploid_metadata, copy=False)
        deme_sizes = np.unique(md['deme'], return_counts=True)
        for i in deme_sizes[1]:
            self.assertEqual(i, self.pop.N//2)
        self.assertTrue(validate_alive_node_metadata(self.pop))

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
            mass_migrations=mmigs, migmatrix=mm, set_selfing_rates=s)
        self.pdict['demography'] = d
        self.pdict['simlen'] = 5
        params = fwdpy11.ModelParams(**self.pdict)
        fwdpy11.evolvets(self.rng, self.pop, params, 100)
        self.assertTrue(validate_alive_node_metadata(self.pop))


class TestGeneticValueLists(unittest.TestCase):
    """
    Added in 0.6.0 to test exceptions raised
    when the number of genetic value objects
    are out of whack vis-a-vis the demography.
    """
    @classmethod
    def setUp(self):
        self.rng = fwdpy11.GSLrng(42)
        self.pop = fwdpy11.DiploidPopulation(100, 1.)
        self.pdict = {'nregions': [],
                      'sregions': [],
                      'recregions': [],
                      'rates': (0, 0, 0),
                      'demography': None,
                      'simlen': 1,
                      'gvalue': fwdpy11.Additive(2.)
                      }

    def test_too_many(self):
        self.pdict['demography'] = np.array([self.pop.N]*10, dtype=np.uint32)
        self.pdict['gvalue'] = [fwdpy11.Additive(2.0)]*2
        params = fwdpy11.ModelParams(**self.pdict)
        with self.assertRaises(ValueError):
            fwdpy11.evolvets(self.rng, self.pop, params, 100)

    def test_too_few(self):
        mmigs = [fwdpy11.move_individuals(
            0, 0, 1, 1./3.), fwdpy11.move_individuals(0, 0, 2, 1./3.)]
        d = fwdpy11.DiscreteDemography(mass_migrations=mmigs)
        self.pdict['demography'] = d
        self.pdict['simlen'] = 5
        self.pdict['gvalue'] = [fwdpy11.Additive(2.0)]*2  # This is the error
        params = fwdpy11.ModelParams(**self.pdict)
        with self.assertRaises(ValueError):
            fwdpy11.evolvets(self.rng, self.pop, params, 100)

    def test_too_few_alt_method(self):
        md = np.array(self.pop.diploid_metadata, copy=False)
        md['deme'][self.pop.N//3:] += 1
        md['deme'][2*self.pop.N//3:] += 1
        self.pdict['demography'] = np.array([self.pop.N]*100, dtype=np.uint32)
        self.pdict['gvalue'] = [fwdpy11.Additive(2.0)]*2  # This is the error
        params = fwdpy11.ModelParams(**self.pdict)
        with self.assertRaises(ValueError):
            fwdpy11.evolvets(self.rng, self.pop, params, 100)


if __name__ == "__main__":
    unittest.main()
