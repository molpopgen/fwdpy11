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

import numpy as np

import fwdpy11


class TestMoveOrCopyIndividuals(unittest.TestCase):
    """
    Tests functions for object construction and
    object pickling
    """

    def test_move(self):
        m = fwdpy11.move_individuals(0, 2, 1, 0.1)
        self.assertEqual(m.when, 0)
        self.assertEqual(m.source, 2)
        self.assertEqual(m.destination, 1)
        self.assertEqual(m.fraction, 0.1)
        self.assertEqual(m.move_individuals, True)
        self.assertEqual(m.resets_growth_rate, True)

        mp = pickle.dumps(m, -1)
        ump = pickle.loads(mp)
        self.assertEqual(ump.when, 0)
        self.assertEqual(ump.source, 2)
        self.assertEqual(ump.destination, 1)
        self.assertEqual(ump.fraction, 0.1)
        self.assertEqual(ump.move_individuals, True)
        self.assertEqual(ump.resets_growth_rate, True)
        self.assertTrue(m == ump)

    def test_copy(self):
        m = fwdpy11.copy_individuals(0, 2, 1, 0.1)
        self.assertEqual(m.when, 0)
        self.assertEqual(m.source, 2)
        self.assertEqual(m.destination, 1)
        self.assertEqual(m.fraction, 0.1)
        self.assertEqual(m.move_individuals, False)

        mp = pickle.dumps(m, -1)
        ump = pickle.loads(mp)
        self.assertEqual(ump.when, 0)
        self.assertEqual(ump.source, 2)
        self.assertEqual(ump.destination, 1)
        self.assertEqual(ump.fraction, 0.1)
        self.assertEqual(ump.move_individuals, False)
        self.assertTrue(m == ump)


class TestSetDemeSize(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.ddsc = fwdpy11.SetDemeSize(0, 1, 3)

    def test_property_access(self):
        self.assertEqual(self.ddsc.when, 0)
        self.assertEqual(self.ddsc.deme, 1)
        self.assertEqual(self.ddsc.new_size, 3)
        self.assertEqual(self.ddsc.resets_growth_rate, True)

    def test_pickle(self):
        p = pickle.dumps(self.ddsc, -1)
        up = pickle.loads(p)
        self.assertEqual(self.ddsc.when, up.when)
        self.assertEqual(self.ddsc.deme, up.deme)
        self.assertEqual(self.ddsc.new_size, up.new_size)

    def test_invalid_deme(self):
        with self.assertRaises(ValueError):
            self.dpsc = fwdpy11.SetDemeSize(0, -1, 3)


class TestSetExponentialGrowth(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.seg = fwdpy11.SetExponentialGrowth(0, 1, 0.1)

    def test_property_access(self):
        self.assertEqual(self.seg.when, 0)
        self.assertEqual(self.seg.deme, 1)
        self.assertEqual(self.seg.G, 0.1)

    def test_pickle(self):
        p = pickle.dumps(self.seg, -1)
        up = pickle.loads(p)
        self.assertEqual(self.seg.when, up.when)
        self.assertEqual(self.seg.deme, up.deme)
        self.assertEqual(self.seg.G, up.G)

    def test_invalid_deme(self):
        with self.assertRaises(ValueError):
            self.seg = fwdpy11.SetExponentialGrowth(0, -1, 0.1)

    def test_invalid_growth_rate(self):
        with self.assertRaises(ValueError):
            self.seg = fwdpy11.SetExponentialGrowth(0, 1, np.nan)


class TestMigrationMatrix(unittest.TestCase):
    def testSimpleInit(self):
        fwdpy11.MigrationMatrix(np.identity(2))

    def testShape(self):
        m = fwdpy11.MigrationMatrix(np.identity(2))
        self.assertEqual(m.shape, (2, 2))

    def testPickle(self):
        m = fwdpy11.MigrationMatrix(np.identity(2))
        p = pickle.dumps(m, -1)
        up = pickle.loads(p)
        self.assertEqual(m, up)
        self.assertEqual(up.shape, (2, 2))
        self.assertTrue(np.array_equal(up.M, np.identity(2)))

    @unittest.skip("obsolete, and a bad idea anyways")
    def testSettingRates(self):
        m = fwdpy11.MigrationMatrix(np.identity(2))
        m._set_migration_rates(0, [0.5, 0.5])
        m._set_migration_rates(1, np.array([0.01, 0.99]))
        self.assertTrue(
            np.array_equal(m.M, np.array([0.5, 0.5, 0.01, 0.99]).reshape(2, 2))
        )

        m._set_migration_rates(np.array([0, 1, 2, 3]).reshape(2, 2))

        try:
            m._set_migration_rates(0, [2.5, 0.5])
        except:  # NOQA
            self.fail("unexpected exception")

        with self.assertRaises(ValueError):
            # Fails because deme index out of range
            m._set_migration_rates(2, [0.5, 0.5])

    def test_non_square_input(self):
        with self.assertRaises(ValueError):
            fwdpy11.MigrationMatrix(np.array([2.0] * 2))

        with self.assertRaises(ValueError):
            """This matrix is also invalid b/c
            rows don't sum to 1, but we expect
            a failure before that check"""
            fwdpy11.MigrationMatrix(np.array([i for i in range(6)]).reshape(2, 3))

    def test_weights_greater_than_one(self):
        with self.assertRaises(ValueError):
            fwdpy11.MigrationMatrix(np.array([2.0] * 4).reshape(2, 2), False)


class TestSetMigrationRates(unittest.TestCase):
    def test_init_from_list(self):
        m = fwdpy11.SetMigrationRates(0, 0, [0, 1, 0])
        self.assertEqual(m.when, 0)
        self.assertEqual(m.deme, 0)
        self.assertTrue(np.array_equal(m.migrates, np.array([0, 1, 0])))

    def test_init_from_numpy(self):
        m = fwdpy11.SetMigrationRates(0, 0, np.array([0, 1, 0]))
        self.assertEqual(m.when, 0)
        self.assertEqual(m.deme, 0)
        self.assertTrue(np.array_equal(m.migrates, np.array([0, 1, 0])))

    def test_reset_entire_matrix(self):
        m = fwdpy11.SetMigrationRates(3, None, np.identity(3))
        self.assertTrue(np.array_equal(m.migrates, np.identity(3)))

    def test_reset_entire_matrix_bad_inputs(self):
        with self.assertRaises(ValueError):
            fwdpy11.SetMigrationRates(3, None, np.arange(4).reshape(2, 2))

    def test_pickle(self):
        m = fwdpy11.SetMigrationRates(0, 0, np.array([0, 1, 0]))
        p = pickle.dumps(m, -1)
        up = pickle.loads(p)
        self.assertEqual(m.when, up.when)
        self.assertEqual(m.deme, up.deme)
        self.assertEqual(m.migrates.tolist(), up.migrates.tolist())


class TestDiscreteDemographyInitialization(unittest.TestCase):
    @unittest.skip("obsolete: UI may need rethinking")
    def test_init_migmatrix_with_tuple(self):
        mm = np.array([0.3, 0.7, 0.7, 0.3]).reshape(2, 2)
        try:
            d = fwdpy11.DiscreteDemography(migmatrix=mm)
            self.assertEqual(d.migmatrix.scaled, False)
        except:  # NOQA
            self.fail("unexpected exception")
        try:
            d = fwdpy11.DiscreteDemography(migmatrix=((mm, True)))
            self.assertEqual(d.migmatrix.scaled, True)
        except:  # NOQA
            self.fail("unexpected exception")
        try:
            d = fwdpy11.DiscreteDemography(migmatrix=((mm, False)))
            self.assertEqual(d.migmatrix.scaled, False)
        except:  # NOQA
            self.fail("unexpected exception")

    def test_move_or_copy(self):
        m = [fwdpy11.move_individuals(0, 0, 1, 0.5)]
        dd = fwdpy11.DiscreteDemography(m)
        self.assertTrue(dd.migmatrix is None)
        with self.assertRaises(AttributeError):
            del self.migmatrix

    def test_setting_growth(self):
        g = [fwdpy11.SetExponentialGrowth(0, 1, 0.3)]
        fwdpy11.DiscreteDemography(set_growth_rates=g)

    def test_setting_size_changes(self):
        c = [fwdpy11.SetDemeSize(0, 1, 15151)]
        fwdpy11.DiscreteDemography(set_deme_sizes=c)

    def test_init_from_numpy(self):
        N = 1000
        popsizes = np.array(
            [N] * 10 * N + [int(0.1 * N)] * int(0.1 * N) + [2 * N] * 100,
            dtype=np.uint32,
        )
        demography = fwdpy11.DiscreteDemography(set_deme_sizes=popsizes)
        self.assertEqual(len(np.unique(popsizes)), len(demography.set_deme_sizes))

    def test_init_with_two_growth_changes_same_generation(self):
        g = [fwdpy11.SetExponentialGrowth(0, 1, 0.3)] * 2
        with self.assertRaises(ValueError):
            fwdpy11.DiscreteDemography(set_growth_rates=g)

    def test_init_with_two_selfing_changes_same_generation(self):
        g = [fwdpy11.SetSelfingRate(0, 1, 0.3)] * 2
        with self.assertRaises(ValueError):
            fwdpy11.DiscreteDemography(set_selfing_rates=g)

    def test_init_with_two_deme_size_changes_same_generation(self):
        c = [fwdpy11.SetDemeSize(0, 1, 15151)] * 2
        with self.assertRaises(ValueError):
            fwdpy11.DiscreteDemography(set_deme_sizes=c)

    def test_init_with_two_mass_moves_same_generation(self):
        m = [fwdpy11.move_individuals(0, 0, 1, 0.5)] * 2
        with self.assertRaises(ValueError):
            fwdpy11.DiscreteDemography(m)

    def test_init_with_two_mass_copies_same_generation(self):
        m = [fwdpy11.copy_individuals(0, 0, 1, 0.5)] * 2
        with self.assertRaises(ValueError):
            fwdpy11.DiscreteDemography(m)

    def test_init_with_one_mass_move_and_one_mass_copy_same_generation(self):
        # Internally, the input data will get sorted
        # so that copies happen before moves.
        m = [
            fwdpy11.move_individuals(0, 0, 1, 0.5),
            fwdpy11.copy_individuals(0, 0, 1, 0.5),
        ]
        try:
            d = fwdpy11.DiscreteDemography(m)
            self.assertEqual(d.mass_migrations[0].move_individuals, False)
            self.assertEqual(d.mass_migrations[1].move_individuals, True)
        except:  # NOQA
            self.fail("unexpected exception")

    def test_move_too_many_individuals_from_same_deme(self):
        # Attempt to move 125% of deme 0
        m = [
            fwdpy11.move_individuals(0, 0, 1, 0.5),
            fwdpy11.move_individuals(0, 0, 2, 0.75),
        ]
        with self.assertRaises(ValueError):
            fwdpy11.DiscreteDemography(m)

    def test_move_too_many_individuals_from_same_deme_with_copy_in_middle(self):
        # Attempt to move 125% of deme 0
        # This is a test that input data get sorted so
        # that copies are before moves, which will allow
        # detecting the invalid cumulative move attempt
        m = [
            fwdpy11.move_individuals(0, 0, 1, 0.5),
            fwdpy11.copy_individuals(0, 0, 1, 0.5),
            fwdpy11.move_individuals(0, 0, 2, 0.75),
        ]
        with self.assertRaises(ValueError):
            fwdpy11.DiscreteDemography(m)

    def test_with_migration_matrix(self):
        mm = np.array([0.72, 0.28, 0.3, 0.7]).reshape(2, 2)
        d = fwdpy11.DiscreteDemography(migmatrix=mm)
        p = pickle.dumps(d, -1)
        up = pickle.loads(p)
        self.assertTrue(np.array_equal(up.migmatrix.M, mm))
        with self.assertRaises(AttributeError):
            del d.migmatrix

    @unittest.skip("obsolete")
    def test_identity_matrix_conversion_to_None(self):
        mm = np.identity(5)
        d = fwdpy11.DiscreteDemography(migmatrix=mm)
        self.assertTrue(d.migmatrix is None)

    def test_identity_matrix_with_migrate_changes(self):
        mm = np.identity(5)
        cm = [fwdpy11.SetMigrationRates(0, 3, [0.2] * 5)]
        d = fwdpy11.DiscreteDemography(migmatrix=mm, set_migration_rates=cm)
        self.assertTrue(d.migmatrix is not None)


class TestEventSorting(unittest.TestCase):
    """
    Doesn't test all possibilities, esp. for MassMigration events
    """

    def test_events_are_sorted(self):
        d = fwdpy11.DiscreteDemography(
            mass_migrations=[
                fwdpy11.MassMigration(when=10, source=0, destination=1, fraction=0.25),
                fwdpy11.MassMigration(when=1, source=1, destination=0, fraction=0.25),
            ],
            set_deme_sizes=[
                fwdpy11.SetDemeSize(when=10, deme=0, new_size=11),
                fwdpy11.SetDemeSize(when=3, deme=1, new_size=100),
            ],
            set_growth_rates=[
                fwdpy11.SetExponentialGrowth(when=10, deme=0, G=0.25),
                fwdpy11.SetExponentialGrowth(when=1, deme=0, G=fwdpy11.NOGROWTH),
            ],
            set_selfing_rates=[
                fwdpy11.SetSelfingRate(when=303, deme=0, S=0.1),
                fwdpy11.SetSelfingRate(when=1, deme=0, S=1.0),
            ],
            set_migration_rates=[
                fwdpy11.SetMigrationRates(when=10, deme=0, migrates=[0.5, 0.5]),
                fwdpy11.SetMigrationRates(when=8, deme=0, migrates=[1, 0.0]),
            ],
            migmatrix=np.eye(2),
        )
        for i in [
            d.mass_migrations,
            d.set_deme_sizes,
            d.set_growth_rates,
            d.set_selfing_rates,
            d.set_migration_rates,
        ]:
            self.assertTrue(i == sorted(i, key=lambda x: x.when))


class TestPickling(unittest.TestCase):
    """
    Test pickling a couple of models that use all of the features.
    """

    def test_pickle_IM_with_selfing(self):
        import fwdpy11.demographic_models.IM as IM

        temp = IM.two_deme_IM(1000, 100, 0.25, [0.5, 3], [0.1, 0.25]).model.asdict()
        temp["set_selfing_rates"] = [fwdpy11.SetSelfingRate(0, 0, 0.1)]
        im = fwdpy11.DiscreteDemography(**temp)
        pi = pickle.dumps(im)
        up = pickle.loads(pi)
        self.assertEqual(im, up)

    def test_pickle_IM_with_whole_migmatrix_reset(self):
        import fwdpy11.demographic_models.IM as IM

        temp = IM.two_deme_IM(1000, 100, 0.25, [0.5, 3], [0.1, 0.25]).model.asdict()
        temp["set_migration_rates"].append(
            fwdpy11.SetMigrationRates(10000 + 10, None, np.eye(2))
        )
        im = fwdpy11.DiscreteDemography(**temp)
        pi = pickle.dumps(im)
        up = pickle.loads(pi)
        self.assertEqual(im, up)

    def test_pickle_tennessen(self):
        import fwdpy11.demographic_models.human as human

        tennessen = human.tennessen()
        pt = pickle.dumps(tennessen)
        up = pickle.loads(pt)
        self.assertEqual(tennessen, up)

    def test_pickle_model_with_scaled_migration_matrix(self):
        m = fwdpy11.DiscreteDemography(migmatrix=(np.eye(2), True))
        p = pickle.dumps(m)
        up = pickle.loads(p)
        self.assertEqual(m, up)


if __name__ == "__main__":
    unittest.main()
