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

import numpy as np

import fwdpy11


def setup_and_run_model(pop, ddemog, simlen, recorder=None, seed=654321):
    pdict = {
        "nregions": [],
        "sregions": [],
        "recregions": [],
        "rates": (0, 0, 0,),
        "gvalue": fwdpy11.Multiplicative(2.0),
        "demography": ddemog,
        "simlen": simlen,
    }
    params = fwdpy11.ModelParams(**pdict)
    rng = fwdpy11.GSLrng(seed)
    fwdpy11.evolvets(rng, pop, params, 100, recorder)


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
        with self.assertRaises(RuntimeError):
            setup_and_run_model(self.pop, d, 5)

    def test_copy_from_empty_deme(self):
        m = [fwdpy11.move_individuals(0, 0, 1, 1), fwdpy11.copy_individuals(1, 0, 1, 1)]
        d = fwdpy11.DiscreteDemography(mass_migrations=m)
        with self.assertRaises(ValueError):
            fwdpy11.DemographyDebugger(self.pop, d)
        with self.assertRaises(RuntimeError):
            setup_and_run_model(self.pop, d, 5)

    def test_failure_to_update_migration_matrix(self):
        M = np.zeros(4).reshape(2, 2)
        M[0, 0] = 1.0
        mm = [fwdpy11.move_individuals(when=10, source=0, destination=1, fraction=0.25)]
        d = fwdpy11.DiscreteDemography(migmatrix=M, mass_migrations=mm)
        with self.assertRaises(ValueError):
            fwdpy11.DemographyDebugger(self.pop, d)


class TestDetectingExtinctions(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.pop = fwdpy11.DiploidPopulation(100, 1.0)

    def test_no_valid_parents_V1(self):
        m = [fwdpy11.move_individuals(0, 0, 1, 1)]
        s = [fwdpy11.SetDemeSize(0, 0, 100)]
        d = fwdpy11.DiscreteDemography(mass_migrations=m, set_deme_sizes=s)
        with self.assertRaises(ValueError):
            fwdpy11.DemographyDebugger(self.pop, d)

    def test_no_valid_parents_V2(self):
        m = [fwdpy11.move_individuals(0, 0, 1, 1)]
        s = [fwdpy11.SetDemeSize(1, 0, 100)]
        d = fwdpy11.DiscreteDemography(mass_migrations=m, set_deme_sizes=s)
        with self.assertRaises(ValueError):
            fwdpy11.DemographyDebugger(self.pop, d)

    def test_no_valid_parents_V3(self):
        m = [fwdpy11.copy_individuals(0, 0, 1, 1)]
        s = [fwdpy11.SetDemeSize(0, 0, 0), fwdpy11.SetDemeSize(11, 0, 100)]
        d = fwdpy11.DiscreteDemography(mass_migrations=m, set_deme_sizes=s)
        with self.assertRaises(ValueError):
            fwdpy11.DemographyDebugger(self.pop, d)

    def test_no_valid_parents_with_migration(self):
        m = [fwdpy11.move_individuals(0, 0, 1, 1)]
        s = [fwdpy11.SetDemeSize(0, 0, 100)]
        M = np.array([0.5] * 4).reshape(2, 2)
        d = fwdpy11.DiscreteDemography(mass_migrations=m, set_deme_sizes=s, migmatrix=M)
        with self.assertRaises(ValueError):
            fwdpy11.DemographyDebugger(self.pop, d)

    def test_no_valid_parents_with_migration_V2(self):
        m = [fwdpy11.move_individuals(0, 0, 1, 1)]
        s = [fwdpy11.SetDemeSize(0, 0, 100)]
        M = np.array([0.5] * 4).reshape(2, 2)
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
        except Exception:
            self.fail("unexpected exception")
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
        except Exception:
            self.fail("unexpected exception")
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
