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

"""
This is basically a test of bad models.
"""

import unittest
import warnings

import numpy as np

import fwdpy11


def setup_and_run_model(pop, ddemog, simlen, recorder=None, seed=654321):
    pdict = {
        "nregions": [],
        "sregions": [],
        "recregions": [],
        "rates": (
            0,
            0,
            0,
        ),
        "gvalue": fwdpy11.Multiplicative(2.0),
        "demography": ddemog,
        "simlen": simlen,
    }
    params = fwdpy11.ModelParams(**pdict)
    rng = fwdpy11.GSLrng(seed)
    fwdpy11.evolvets(rng, pop, params, 100, recorder)


class TestNoDemographicEvents(unittest.TestCase):
    """
    For issue #594
    """

    @classmethod
    def setUpClass(self):
        self.demog = fwdpy11.DiscreteDemography()

    def test_no_demographic_events(self):
        _ = fwdpy11.DemographyDebugger([100], self.demog)


class TestEmptyEventLists(unittest.TestCase):
    """
    For issue #594
    """

    @classmethod
    def setUpClass(self):
        self.demog = fwdpy11.DiscreteDemography(
            mass_migrations=[],
            set_growth_rates=[],
            set_deme_sizes=[],
            set_selfing_rates=[],
            migmatrix=None,
            set_migration_rates=[],
        )

    def test_empty_event_lists(self):
        _ = fwdpy11.DemographyDebugger([100], self.demog)


class TestSingleDemeWithVariousMigrationMatrices(unittest.TestCase):
    """
    For issue #594
    """

    @classmethod
    def setUp(self):
        self.pop = fwdpy11.DiploidPopulation(100, 1.0)

    def test_one_by_one_migration_matrix(self):
        """
        Will warn about the 1-by-1 matrix
        """
        demog = fwdpy11.DiscreteDemography(migmatrix=np.array([[1.0]]))
        warned = False
        with warnings.catch_warnings(record=True) as w:
            _ = fwdpy11.DemographyDebugger([100], demog)
            warned = True
        assert warned is True

    def test_two_by_two_migration_matrix(self):
        demog = fwdpy11.DiscreteDemography(migmatrix=np.array([[1.0, 0.0], [0.0, 1.0]]))
        with self.assertRaises(ValueError):
            _ = fwdpy11.DemographyDebugger([100], demog)


class TestBadMassMigrations(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.pop = fwdpy11.DiploidPopulation(100, 1.0)

    def test_move_too_many(self):
        m = [
            fwdpy11.move_individuals(0, 0, 1, 0.5),
            fwdpy11.move_individuals(1, 0, 1, 1),
            fwdpy11.move_individuals(2, 0, 1, 1),
        ]
        d = fwdpy11.DiscreteDemography(mass_migrations=m)
        with self.assertRaises(ValueError):
            fwdpy11.DemographyDebugger(self.pop, d)
        with self.assertRaises(fwdpy11.DemographyError):
            setup_and_run_model(self.pop, d, 5)

    def test_copy_from_empty_deme(self):
        m = [fwdpy11.move_individuals(0, 0, 1, 1), fwdpy11.copy_individuals(1, 0, 1, 1)]
        d = fwdpy11.DiscreteDemography(mass_migrations=m)
        with self.assertRaises(ValueError):
            fwdpy11.DemographyDebugger(self.pop, d)
        with self.assertRaises(fwdpy11.DemographyError):
            setup_and_run_model(self.pop, d, 5)

    def test_failure_to_update_migration_matrix(self):
        M = np.zeros(4).reshape(2, 2)
        M[0, 0] = 1.0
        mm = [fwdpy11.move_individuals(when=10, source=0, destination=1, fraction=0.25)]
        d = fwdpy11.DiscreteDemography(migmatrix=M, mass_migrations=mm)
        with self.assertRaises(ValueError):
            fwdpy11.DemographyDebugger(self.pop, d)
        with self.assertRaises(fwdpy11.DemographyError):
            setup_and_run_model(self.pop, d, 20)


class TestDetectingAbsenceOfAncestry(unittest.TestCase):
    """
    These tests check for errors arising due to
    there being no valid parental generation at some
    point.
    """
    @classmethod
    def setUp(self):
        self.pop = fwdpy11.DiploidPopulation(100, 1.0)

    def test_no_valid_parents_V1(self):
        m = [fwdpy11.move_individuals(0, 0, 1, 1)]
        s = [fwdpy11.SetDemeSize(0, 0, 100)]
        d = fwdpy11.DiscreteDemography(mass_migrations=m, set_deme_sizes=s)
        with self.assertRaises(ValueError):
            fwdpy11.DemographyDebugger(self.pop, d)
        with self.assertRaises(fwdpy11.DemographyError):
            setup_and_run_model(self.pop, d, 5)

    def test_no_valid_parents_V2(self):
        m = [fwdpy11.move_individuals(0, 0, 1, 1)]
        s = [fwdpy11.SetDemeSize(1, 0, 100)]
        d = fwdpy11.DiscreteDemography(mass_migrations=m, set_deme_sizes=s)
        with self.assertRaises(ValueError):
            fwdpy11.DemographyDebugger(self.pop, d)
        with self.assertRaises(fwdpy11.DemographyError):
            setup_and_run_model(self.pop, d, 5)

    def test_no_valid_parents_V3(self):
        m = [fwdpy11.copy_individuals(0, 0, 1, 1)]
        s = [fwdpy11.SetDemeSize(0, 0, 0), fwdpy11.SetDemeSize(11, 0, 100)]
        d = fwdpy11.DiscreteDemography(mass_migrations=m, set_deme_sizes=s)
        with self.assertRaises(ValueError):
            fwdpy11.DemographyDebugger(self.pop, d)
        with self.assertRaises(fwdpy11.DemographyError):
            setup_and_run_model(self.pop, d, 12)

    def test_no_valid_parents_with_migration(self):
        m = [fwdpy11.move_individuals(0, 0, 1, 1)]
        s = [fwdpy11.SetDemeSize(0, 0, 100)]
        M = np.array([0.5] * 4).reshape(2, 2)
        d = fwdpy11.DiscreteDemography(mass_migrations=m, set_deme_sizes=s, migmatrix=M)
        with self.assertRaises(ValueError):
            fwdpy11.DemographyDebugger(self.pop, d)
        with self.assertRaises(fwdpy11.DemographyError):
            setup_and_run_model(self.pop, d, 5)

    def test_no_valid_parents_with_migration_V2(self):
        """
        This is a complex test: we mess around with
        things at time 0, which is relevent to scenarios
        where we stop/start models.

        A complementary tests exists in the C++ suite.
        """
        # Move all founders from deme 0 to deme 1
        # In gen 1, all ancestry in deme 1 is from
        # deme 0
        m = [fwdpy11.move_individuals(0, 0, 1, 1)]
        # At the same time, set deme 0 size to 100
        s = [fwdpy11.SetDemeSize(0, 0, 100)]
        # The initial migration matrix will have 50%
        # ancestry for all pops each generation...
        M = np.array([0.5] * 4).reshape(2, 2)
        # But we will immediately change this
        # so that all ancestry for both demes
        # is from deme 1, making the model valid
        cM = [
            fwdpy11.SetMigrationRates(
                when=0, deme=None, migrates=np.array([0, 1, 0, 1]).reshape(2, 2)
            )
        ]
        d = fwdpy11.DiscreteDemography(
            mass_migrations=m, set_deme_sizes=s, set_migration_rates=cM, migmatrix=M
        )
        try:
            fwdpy11.DemographyDebugger(self.pop, d)
        except Exception as e:
            self.fail(f"unexpected exception, {e}")
        setup_and_run_model(self.pop, d, 5)

    def test_no_valid_parents_with_migration_V2_migration_rescue(self):
        """
        This is a convoluted model:
        0 moves 100% to 1. Thus, there are no parents anymore for 0.
        But, this is okay because the new deme 1 will provide parents
        back into 0, thus "un-extincting" 0.  Then, in the next
        generation, each deme gets 50/50 of its ancestry from each
        deme.
        """
        m = [fwdpy11.move_individuals(0, 0, 1, 1)]
        s = [fwdpy11.SetDemeSize(0, 0, 100)]
        M = np.array([0.5] * 4).reshape(2, 2)
        cM = [
            fwdpy11.SetMigrationRates(0, None, np.array([0, 1, 0, 1]).reshape(2, 2)),
            fwdpy11.SetMigrationRates(1, None, np.array([0.5] * 4).reshape(2, 2)),
        ]
        d = fwdpy11.DiscreteDemography(
            mass_migrations=m, set_deme_sizes=s, set_migration_rates=cM, migmatrix=M
        )
        try:
            fwdpy11.DemographyDebugger(self.pop, d)
        except Exception as e:
            self.fail(f"unexpected exception, {e}")
        setup_and_run_model(self.pop, d, 5)

    def test_no_valid_parents_with_migration_V3(self):
        m = [fwdpy11.move_individuals(0, 0, 1, 1)]
        s = [fwdpy11.SetDemeSize(0, 0, 100)]
        M = np.array([0.5] * 4).reshape(2, 2)
        d = fwdpy11.DiscreteDemography(
            mass_migrations=m,
            set_deme_sizes=s,
            migmatrix=fwdpy11.MigrationMatrix(M, True),
        )
        with self.assertRaises(ValueError):
            fwdpy11.DemographyDebugger(self.pop, d)
        with self.assertRaises(fwdpy11.DemographyError):
            setup_and_run_model(self.pop, d, 5)


class TestBadModels(unittest.TestCase):
    def test_github_issue_544(self):
        pop = fwdpy11.DiploidPopulation(100, 1)
        mass_migrations = [
            fwdpy11.move_individuals(when=100, source=0, destination=1, fraction=0.5)
        ]
        m = 1e-3
        migmatrix = np.array([1.0 - m, m, m, 1.0 - m]).reshape(2, 2)
        bad_dmodel = fwdpy11.DiscreteDemography(
            mass_migrations=mass_migrations, migmatrix=migmatrix
        )
        with self.assertRaises(ValueError):
            fwdpy11.DemographyDebugger([100], bad_dmodel)

        pdict = {
            "nregions": [],
            "sregions": [],
            "recregions": [],
            "gvalue": fwdpy11.Multiplicative(1.0),
            "rates": [0, 0, 0],
            "demography": bad_dmodel,
            "simlen": 150,
        }

        rng = fwdpy11.GSLrng(1010)
        mp = fwdpy11.ModelParams(**pdict)
        with self.assertRaises(fwdpy11.DemographyError):
            fwdpy11.evolvets(rng, pop, mp, 100)


if __name__ == "__main__":
    unittest.main()
