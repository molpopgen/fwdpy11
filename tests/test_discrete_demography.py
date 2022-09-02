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
import typing
import glob
import unittest
from dataclasses import dataclass

import attr
import fwdpy11
import intervaltree
import numpy as np
import pytest


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


@pytest.fixture(params=[10234, 91246412, 8513251, 51285123, 252101511])
def random_migrates_seed(request):
    return request.param


@pytest.fixture(params=[i for i in range(2, 30, 2)])
def random_migrates_npops(request):
    return request.param


def test_set_migration_rates_numerical_tolerance(
    random_migrates_seed, random_migrates_npops
):
    """
    Test related to GitHub issue 787
    """
    nreps = 100
    rng = np.random.Generator(np.random.MT19937(seed=random_migrates_seed))
    for index in range(random_migrates_npops):
        alpha = [rng.integers(1, high=10)] * random_migrates_npops
        samples = rng.dirichlet(alpha, size=nreps)
        for row in range(nreps):
            r = samples[row, :][:]
            r[index] = 0.0
            r[index] = 1.0 - r.sum()
            _ = fwdpy11.SetMigrationRates(when=1, deme=index, migrates=r)


# NOTE: from here down we have tests of non-public stuff


# NOTE: use attr b/c we want kw_only, which isn't in dataclass until 3.10
@attr.s(kw_only=True, auto_attribs=True)
class Payload:
    initial_sizes: typing.Dict[int, int]
    mass_migrations: typing.Optional[typing.List[fwdpy11.MassMigration]]
    set_deme_sizes: typing.Optional[typing.List[fwdpy11.SetDemeSize]]
    set_growth_rates: typing.Optional[typing.List[fwdpy11.SetExponentialGrowth]]
    total_simulation_length: typing.Optional[int]
    expected_itree: intervaltree.IntervalTree


def __detailed_size_history_test_runner(payload: Payload):
    import fwdpy11.discrete_demography

    deme_properties = fwdpy11.discrete_demography._DemeSizeHistory.from_lowlevel(
        initial_sizes=payload.initial_sizes,
        mass_migrations=payload.mass_migrations,
        set_deme_sizes=payload.set_deme_sizes,
        set_growth_rates=payload.set_growth_rates,
        total_simulation_length=payload.total_simulation_length,
    )

    assert deme_properties.epochs == payload.expected_itree

    @dataclass
    class TrackSizes:
        sizes: typing.List[typing.Tuple[int, typing.Dict[int, int]]]

        def __call__(self, pop, _):
            self.sizes.append((pop.generation, pop.deme_sizes(as_dict=True)))

    t = TrackSizes([])

    demog = fwdpy11.DiscreteDemography(
        set_deme_sizes=payload.set_deme_sizes,
        mass_migrations=payload.mass_migrations,
        set_growth_rates=payload.set_growth_rates,
    )

    pdict = {
        "rates": (0, 0, 0),
        "gvalue": fwdpy11.Multiplicative(2.0),
        "simlen": payload.total_simulation_length,
        "demography": demog,
    }

    if pdict["simlen"] is None:
        pdict["simlen"] = 50

    params = fwdpy11.ModelParams(**pdict)

    pop = fwdpy11.DiploidPopulation([v for _, v in payload.initial_sizes.items()], 1.0)

    rng = fwdpy11.GSLrng(42)

    final_time = params.simlen

    try:
        fwdpy11.evolvets(rng, pop, params, 100, t)
    except fwdpy11.GlobalExtinction as _:
        final_time = pop.generation
    except Exception as e:
        pytest.fail(f"unexpected exception, {e}")

    # The data from the forward sim must overlap an interval
    # in the interval tree
    for i in t.sizes:
        generation, sizes = i
        for deme, size in sizes.items():
            assert deme_properties.deme_exists_at(
                deme, generation
            ), f"failure at {generation} for deme {deme}, which has size {size}"
            # Get the interval for this deme
            for i in deme_properties.epochs[generation]:
                if i.data.deme == deme and i.begin == generation:
                    assert size == int(
                        np.rint(i.data.initial_size * i.data.growth_parameter)
                    )
            predicted_size = deme_properties.deme_size_at(deme, generation)
            assert size is not None
            assert (
                size == predicted_size
            ), f"{size} != {predicted_size} at time {generation}"

    # Check the other direction--does the interval tree predict the sim outputs?
    generations = [0] * (final_time + 1)
    for i in deme_properties.epochs:
        if i.begin == 0:
            start = 1
        else:
            start = i.begin
        for g in range(start, min(i.end, params.simlen + 1)):
            generations[g] += 1
            record = [j[1] for j in t.sizes if j[0] == g]
            assert len(record) == 1
            assert i.data.deme in [key for key in record[0].keys()], f"{record[0]}"
    assert all([i > 0 for i in generations[1:]])


@pytest.mark.parametrize(
    "payload",
    [
        # Relatively simple models
        Payload(
            initial_sizes={0: 10},
            mass_migrations=None,
            set_deme_sizes=None,
            set_growth_rates=None,
            total_simulation_length=None,
            expected_itree=intervaltree.IntervalTree(
                [
                    intervaltree.Interval(
                        0,
                        int(np.iinfo(np.uint32).max),
                        fwdpy11.discrete_demography._EpochData(0, [0], 10, 1.0),
                    )
                ]
            ),
        ),
        Payload(
            initial_sizes={0: 10, 1: 100},
            mass_migrations=None,
            set_deme_sizes=[fwdpy11.SetDemeSize(when=11, deme=0, new_size=777)],
            set_growth_rates=None,
            total_simulation_length=500,
            expected_itree=intervaltree.IntervalTree(
                [
                    intervaltree.Interval(
                        0, 12, fwdpy11.discrete_demography._EpochData(0, [0], 10, 1.0)
                    ),
                    intervaltree.Interval(
                        12,
                        501,
                        fwdpy11.discrete_demography._EpochData(0, [0], 777, 1.0),
                    ),
                    intervaltree.Interval(
                        0, 501, fwdpy11.discrete_demography._EpochData(1, [1], 100, 1.0)
                    ),
                ]
            ),
        ),
        Payload(
            initial_sizes={0: 10, 1: 100},
            mass_migrations=None,
            set_deme_sizes=[
                fwdpy11.SetDemeSize(when=11, deme=0, new_size=777),
                fwdpy11.SetDemeSize(when=14, deme=1, new_size=0),
            ],
            set_growth_rates=None,
            total_simulation_length=500,
            expected_itree=intervaltree.IntervalTree(
                [
                    intervaltree.Interval(
                        0, 12, fwdpy11.discrete_demography._EpochData(0, [0], 10, 1.0)
                    ),
                    intervaltree.Interval(
                        12,
                        501,
                        fwdpy11.discrete_demography._EpochData(0, [0], 777, 1.0),
                    ),
                    intervaltree.Interval(
                        0, 15, fwdpy11.discrete_demography._EpochData(1, [1], 100, 1.0)
                    ),
                ]
            ),
        ),
        Payload(
            initial_sizes={0: 10, 1: 100},
            mass_migrations=None,
            set_deme_sizes=[
                fwdpy11.SetDemeSize(when=11, deme=0, new_size=0),
                fwdpy11.SetDemeSize(when=14, deme=1, new_size=0),
            ],
            set_growth_rates=None,
            total_simulation_length=500,
            expected_itree=intervaltree.IntervalTree(
                [
                    intervaltree.Interval(
                        0, 12, fwdpy11.discrete_demography._EpochData(0, [0], 10, 1.0)
                    ),
                    intervaltree.Interval(
                        0, 15, fwdpy11.discrete_demography._EpochData(1, [1], 100, 1.0)
                    ),
                ]
            ),
        ),
        Payload(
            initial_sizes={0: 10, 1: 100},
            mass_migrations=None,
            set_deme_sizes=[
                fwdpy11.SetDemeSize(when=11, deme=0, new_size=0),
                fwdpy11.SetDemeSize(when=14, deme=1, new_size=0),
            ],
            set_growth_rates=None,
            total_simulation_length=None,
            expected_itree=intervaltree.IntervalTree(
                [
                    intervaltree.Interval(
                        0, 12, fwdpy11.discrete_demography._EpochData(0, [0], 10, 1.0)
                    ),
                    intervaltree.Interval(
                        0, 15, fwdpy11.discrete_demography._EpochData(1, [1], 100, 1.0)
                    ),
                ]
            ),
        ),
        # Simple growth in a single deme
        Payload(
            initial_sizes={0: 10},
            mass_migrations=None,
            set_deme_sizes=None,
            set_growth_rates=[fwdpy11.SetExponentialGrowth(deme=0, when=11, G=1.01)],
            total_simulation_length=None,
            expected_itree=intervaltree.IntervalTree(
                [
                    intervaltree.Interval(
                        0,
                        12,
                        fwdpy11.discrete_demography._EpochData(0, [0], 10, 1.0),
                    ),
                    intervaltree.Interval(
                        12,
                        int(np.iinfo(np.uint32).max),
                        fwdpy11.discrete_demography._EpochData(0, [0], 10, 1.01),
                    ),
                ]
            ),
        ),
        Payload(
            initial_sizes={0: 10},
            mass_migrations=None,
            set_deme_sizes=None,
            set_growth_rates=[
                fwdpy11.SetExponentialGrowth(deme=0, when=1, G=1.0),
            ],
            total_simulation_length=15,
            expected_itree=intervaltree.IntervalTree(
                [
                    intervaltree.Interval(
                        0,
                        2,
                        fwdpy11.discrete_demography._EpochData(0, [0], 10, 1.0),
                    ),
                    intervaltree.Interval(
                        2,
                        16,
                        fwdpy11.discrete_demography._EpochData(0, [0], 10, 1.00),
                    ),
                ]
            ),
        ),
        Payload(
            initial_sizes={0: 10},
            mass_migrations=None,
            set_deme_sizes=None,
            set_growth_rates=[
                fwdpy11.SetExponentialGrowth(deme=0, when=0, G=1.1),
                fwdpy11.SetExponentialGrowth(deme=0, when=1, G=1.0),
            ],
            total_simulation_length=15,
            expected_itree=intervaltree.IntervalTree(
                [
                    intervaltree.Interval(
                        0,
                        1,
                        fwdpy11.discrete_demography._EpochData(0, [0], 10, 1.0),
                    ),
                    intervaltree.Interval(
                        1,
                        2,
                        fwdpy11.discrete_demography._EpochData(0, [0], 10, 1.1),
                    ),
                    intervaltree.Interval(
                        2,
                        16,
                        fwdpy11.discrete_demography._EpochData(0, [0], 11, 1.00),
                    ),
                ]
            ),
        ),
        Payload(
            initial_sizes={0: 10},
            mass_migrations=None,
            set_deme_sizes=None,
            set_growth_rates=[
                fwdpy11.SetExponentialGrowth(deme=0, when=5, G=1.1),
                fwdpy11.SetExponentialGrowth(deme=0, when=10, G=1.02),
            ],
            total_simulation_length=15,
            expected_itree=intervaltree.IntervalTree(
                [
                    intervaltree.Interval(
                        0,
                        6,
                        fwdpy11.discrete_demography._EpochData(0, [0], 10, 1.0),
                    ),
                    intervaltree.Interval(
                        6,
                        11,
                        fwdpy11.discrete_demography._EpochData(0, [0], 10, 1.1),
                    ),
                    intervaltree.Interval(
                        11,
                        16,
                        fwdpy11.discrete_demography._EpochData(0, [0], 16, 1.02),
                    ),
                ]
            ),
        ),
        # Same as above, but we do a hard size reset
        Payload(
            initial_sizes={0: 10},
            mass_migrations=None,
            set_deme_sizes=[fwdpy11.SetDemeSize(when=10, deme=0, new_size=20)],
            set_growth_rates=[
                fwdpy11.SetExponentialGrowth(deme=0, when=5, G=1.1),
                fwdpy11.SetExponentialGrowth(deme=0, when=10, G=1.02),
            ],
            total_simulation_length=15,
            expected_itree=intervaltree.IntervalTree(
                [
                    intervaltree.Interval(
                        0,
                        6,
                        fwdpy11.discrete_demography._EpochData(0, [0], 10, 1.0),
                    ),
                    intervaltree.Interval(
                        6,
                        11,
                        fwdpy11.discrete_demography._EpochData(0, [0], 10, 1.1),
                    ),
                    intervaltree.Interval(
                        11,
                        16,
                        fwdpy11.discrete_demography._EpochData(0, [0], 20, 1.02),
                    ),
                ]
            ),
        ),
        # Same as above, but we do a hard size reset that DOES NOT affect growth rates
        Payload(
            initial_sizes={0: 10},
            mass_migrations=None,
            set_deme_sizes=[
                fwdpy11.SetDemeSize(
                    when=10, deme=0, new_size=20, resets_growth_rate=False
                )
            ],
            set_growth_rates=[
                fwdpy11.SetExponentialGrowth(deme=0, when=5, G=1.1),
            ],
            total_simulation_length=15,
            expected_itree=intervaltree.IntervalTree(
                [
                    intervaltree.Interval(
                        0,
                        6,
                        fwdpy11.discrete_demography._EpochData(0, [0], 10, 1.0),
                    ),
                    intervaltree.Interval(
                        6,
                        11,
                        fwdpy11.discrete_demography._EpochData(0, [0], 10, 1.1),
                    ),
                    intervaltree.Interval(
                        11,
                        16,
                        fwdpy11.discrete_demography._EpochData(0, [0], 20, 1.1),
                    ),
                ]
            ),
        ),
        # Quirky edge cases go here (with comments)
        # One deme, goes extinct at the end of generation 0
        Payload(
            initial_sizes={0: 10},
            mass_migrations=None,
            set_deme_sizes=[fwdpy11.SetDemeSize(deme=0, when=1, new_size=0)],
            set_growth_rates=None,
            total_simulation_length=None,
            expected_itree=intervaltree.IntervalTree(
                [
                    intervaltree.Interval(
                        0,
                        2,
                        fwdpy11.discrete_demography._EpochData(0, [0], 10, 1.0),
                    )
                ]
            ),
        ),
    ],
)
def test_deme_size_history(payload: Payload):
    __detailed_size_history_test_runner(payload)


@pytest.mark.parametrize(
    "payload",
    [
        Payload(
            initial_sizes={0: 10},
            mass_migrations=[
                fwdpy11.copy_individuals(when=10, source=0, destination=1, fraction=1.0)
            ],
            set_deme_sizes=None,
            set_growth_rates=None,
            total_simulation_length=20,
            expected_itree=intervaltree.IntervalTree(
                [
                    intervaltree.Interval(
                        0,
                        21,
                        fwdpy11.discrete_demography._EpochData(0, [0], 10, 1.0),
                    ),
                    intervaltree.Interval(
                        11,
                        21,
                        fwdpy11.discrete_demography._EpochData(1, [0], 10, 1.0),
                    ),
                ]
            ),
        ),
        Payload(
            initial_sizes={0: 100},
            mass_migrations=[
                fwdpy11.copy_individuals(
                    when=10, source=0, destination=1, fraction=0.5
                ),
                fwdpy11.copy_individuals(
                    when=20, source=1, destination=0, fraction=1.0
                ),
            ],
            set_deme_sizes=None,
            set_growth_rates=None,
            total_simulation_length=30,
            expected_itree=intervaltree.IntervalTree(
                [
                    intervaltree.Interval(
                        0,
                        21,
                        fwdpy11.discrete_demography._EpochData(0, [0], 100, 1.0),
                    ),
                    intervaltree.Interval(
                        21,
                        31,
                        fwdpy11.discrete_demography._EpochData(0, [0, 1], 150, 1.0),
                    ),
                    intervaltree.Interval(
                        11,
                        31,
                        fwdpy11.discrete_demography._EpochData(1, [0], 50, 1.0),
                    ),
                ]
            ),
        ),
        # Same as above, but with move in lieu of copies
        Payload(
            initial_sizes={0: 100},
            mass_migrations=[
                fwdpy11.move_individuals(
                    when=10, source=0, destination=1, fraction=0.5
                ),
                fwdpy11.move_individuals(
                    when=20, source=1, destination=0, fraction=1.0
                ),
            ],
            set_deme_sizes=None,
            set_growth_rates=None,
            total_simulation_length=30,
            expected_itree=intervaltree.IntervalTree(
                [
                    intervaltree.Interval(
                        0,
                        11,
                        fwdpy11.discrete_demography._EpochData(0, [0], 100, 1.0),
                    ),
                    intervaltree.Interval(
                        11,
                        21,
                        fwdpy11.discrete_demography._EpochData(0, [0], 50, 1.0),
                    ),
                    intervaltree.Interval(
                        21,
                        31,
                        fwdpy11.discrete_demography._EpochData(0, [0, 1], 100, 1.0),
                    ),
                    intervaltree.Interval(
                        11,
                        21,
                        fwdpy11.discrete_demography._EpochData(1, [0], 50, 1.0),
                    ),
                ]
            ),
        ),
        # Mass migrations and growth.
        # This is where U+1F4A9 hits the fan
        Payload(
            initial_sizes={0: 10},
            mass_migrations=[
                fwdpy11.copy_individuals(
                    when=10, source=0, destination=1, fraction=0.5
                ),
                fwdpy11.copy_individuals(
                    when=15, source=1, destination=0, fraction=1.0
                ),
            ],
            set_deme_sizes=None,
            set_growth_rates=[fwdpy11.SetExponentialGrowth(when=10, deme=1, G=2)],
            total_simulation_length=20,
            expected_itree=intervaltree.IntervalTree(
                [
                    intervaltree.Interval(
                        0,
                        16,
                        fwdpy11.discrete_demography._EpochData(0, [0], 10, 1.0),
                    ),
                    intervaltree.Interval(
                        16,
                        21,
                        fwdpy11.discrete_demography._EpochData(0, [0, 1], 170, 1.0),
                    ),
                    intervaltree.Interval(
                        11,
                        16,
                        fwdpy11.discrete_demography._EpochData(1, [0], 5, 2.0),
                    ),
                    intervaltree.Interval(
                        16,
                        21,
                        fwdpy11.discrete_demography._EpochData(1, [0, 1], 160, 1.0),
                    ),
                ]
            ),
        ),
        # Rare use case: mass migrations
        # not affecting growth rates
        Payload(
            initial_sizes={0: 10},
            mass_migrations=[
                fwdpy11.copy_individuals(
                    when=10,
                    source=0,
                    destination=1,
                    fraction=0.5,
                    resets_growth_rate=False,
                ),
                fwdpy11.copy_individuals(
                    when=15,
                    source=1,
                    destination=0,
                    fraction=1.0,
                    resets_growth_rate=False,
                ),
            ],
            set_deme_sizes=None,
            set_growth_rates=[fwdpy11.SetExponentialGrowth(when=10, deme=1, G=2)],
            total_simulation_length=20,
            expected_itree=intervaltree.IntervalTree(
                [
                    intervaltree.Interval(
                        0,
                        16,
                        fwdpy11.discrete_demography._EpochData(0, [0], 10, 1.0),
                    ),
                    intervaltree.Interval(
                        16,
                        21,
                        fwdpy11.discrete_demography._EpochData(0, [0, 1], 170, 1.0),
                    ),
                    intervaltree.Interval(
                        11,
                        16,
                        fwdpy11.discrete_demography._EpochData(1, [0], 5, 2.0),
                    ),
                    intervaltree.Interval(
                        16,
                        21,
                        fwdpy11.discrete_demography._EpochData(1, [0, 1], 160, 2.0),
                    ),
                ]
            ),
        ),
        Payload(
            initial_sizes={0: 10},
            mass_migrations=[
                fwdpy11.copy_individuals(
                    when=10, source=0, destination=1, fraction=0.5
                ),
                fwdpy11.copy_individuals(
                    when=11, source=1, destination=0, fraction=1.0
                ),
            ],
            set_deme_sizes=None,
            set_growth_rates=None,
            total_simulation_length=20,
            expected_itree=intervaltree.IntervalTree(
                [
                    intervaltree.Interval(
                        0,
                        12,
                        fwdpy11.discrete_demography._EpochData(0, [0], 10, 1.0),
                    ),
                    intervaltree.Interval(
                        12,
                        21,
                        fwdpy11.discrete_demography._EpochData(0, [0, 1], 15, 1.0),
                    ),
                    intervaltree.Interval(
                        11,
                        21,
                        fwdpy11.discrete_demography._EpochData(1, [0], 5, 1.0),
                    ),
                ]
            ),
        ),
        Payload(
            initial_sizes={0: 10},
            mass_migrations=[
                fwdpy11.move_individuals(
                    when=10, source=0, destination=1, fraction=0.5
                ),
                fwdpy11.move_individuals(
                    when=11, source=1, destination=0, fraction=1.0
                ),
            ],
            set_deme_sizes=None,
            set_growth_rates=None,
            total_simulation_length=20,
            expected_itree=intervaltree.IntervalTree(
                [
                    intervaltree.Interval(
                        0,
                        11,
                        fwdpy11.discrete_demography._EpochData(0, [0], 10, 1.0),
                    ),
                    intervaltree.Interval(
                        11,
                        12,
                        fwdpy11.discrete_demography._EpochData(0, [0], 5, 1.0),
                    ),
                    intervaltree.Interval(
                        12,
                        21,
                        fwdpy11.discrete_demography._EpochData(0, [0, 1], 10, 1.0),
                    ),
                    intervaltree.Interval(
                        11,
                        12,
                        fwdpy11.discrete_demography._EpochData(1, [0], 5, 1.0),
                    ),
                ]
            ),
        ),
        Payload(
            initial_sizes={0: 10, 1: 10},
            mass_migrations=[
                fwdpy11.copy_individuals(when=1, source=1, destination=0, fraction=1.0),
            ],
            set_deme_sizes=None,
            set_growth_rates=None,
            total_simulation_length=20,
            expected_itree=intervaltree.IntervalTree(
                [
                    intervaltree.Interval(
                        0,
                        2,
                        fwdpy11.discrete_demography._EpochData(0, [0], 10, 1.0),
                    ),
                    intervaltree.Interval(
                        2,
                        21,
                        fwdpy11.discrete_demography._EpochData(0, [0, 1], 20, 1.0),
                    ),
                    intervaltree.Interval(
                        0, 21, fwdpy11.discrete_demography._EpochData(1, [1], 10, 1.0)
                    ),
                ]
            ),
        ),
        Payload(
            initial_sizes={0: 10, 1: 10},
            mass_migrations=[
                fwdpy11.copy_individuals(when=1, source=1, destination=0, fraction=1.0),
            ],
            set_deme_sizes=[fwdpy11.SetDemeSize(when=5, deme=0, new_size=10)],
            set_growth_rates=None,
            total_simulation_length=20,
            expected_itree=intervaltree.IntervalTree(
                [
                    intervaltree.Interval(
                        0,
                        2,
                        fwdpy11.discrete_demography._EpochData(0, [0], 10, 1.0),
                    ),
                    intervaltree.Interval(
                        2,
                        6,
                        fwdpy11.discrete_demography._EpochData(0, [0, 1], 20, 1.0),
                    ),
                    intervaltree.Interval(
                        6,
                        21,
                        fwdpy11.discrete_demography._EpochData(0, [0], 10, 1.0),
                    ),
                    intervaltree.Interval(
                        0, 21, fwdpy11.discrete_demography._EpochData(1, [1], 10, 1.0)
                    ),
                ]
            ),
        ),
    ],
)
def test_deme_size_history_two_deme_mass_migrations(payload: Payload):
    __detailed_size_history_test_runner(payload)


@pytest.mark.parametrize(
    "payload",
    [
        Payload(
            initial_sizes={0: 10, 1: 10, 2: 10},
            mass_migrations=[
                fwdpy11.copy_individuals(when=1, source=1, destination=0, fraction=1.0),
                fwdpy11.copy_individuals(when=1, source=2, destination=0, fraction=1.0),
            ],
            set_deme_sizes=None,
            set_growth_rates=None,
            total_simulation_length=20,
            expected_itree=intervaltree.IntervalTree(
                [
                    intervaltree.Interval(
                        0,
                        2,
                        fwdpy11.discrete_demography._EpochData(0, [0], 10, 1.0),
                    ),
                    intervaltree.Interval(
                        2,
                        21,
                        fwdpy11.discrete_demography._EpochData(0, [0, 1, 2], 30, 1.0),
                    ),
                    intervaltree.Interval(
                        0, 21, fwdpy11.discrete_demography._EpochData(1, [1], 10, 1.0)
                    ),
                    intervaltree.Interval(
                        0, 21, fwdpy11.discrete_demography._EpochData(2, [2], 10, 1.0)
                    ),
                ]
            ),
        ),
        Payload(
            initial_sizes={0: 10, 1: 10, 2: 10},
            mass_migrations=[
                fwdpy11.copy_individuals(when=1, source=1, destination=0, fraction=1.0),
                fwdpy11.copy_individuals(when=1, source=2, destination=0, fraction=1.0),
                fwdpy11.copy_individuals(when=5, source=2, destination=0, fraction=1.0),
            ],
            set_deme_sizes=None,
            set_growth_rates=None,
            total_simulation_length=20,
            expected_itree=intervaltree.IntervalTree(
                [
                    intervaltree.Interval(
                        0,
                        2,
                        fwdpy11.discrete_demography._EpochData(0, [0], 10, 1.0),
                    ),
                    intervaltree.Interval(
                        2,
                        6,
                        fwdpy11.discrete_demography._EpochData(0, [0, 1, 2], 30, 1.0),
                    ),
                    intervaltree.Interval(
                        6,
                        21,
                        fwdpy11.discrete_demography._EpochData(0, [0, 2], 40, 1.0),
                    ),
                    intervaltree.Interval(
                        0, 21, fwdpy11.discrete_demography._EpochData(1, [1], 10, 1.0)
                    ),
                    intervaltree.Interval(
                        0, 21, fwdpy11.discrete_demography._EpochData(2, [2], 10, 1.0)
                    ),
                ]
            ),
        ),
        Payload(
            initial_sizes={0: 10, 1: 10, 2: 10},
            mass_migrations=[
                fwdpy11.copy_individuals(when=1, source=1, destination=0, fraction=1.0),
                fwdpy11.copy_individuals(when=1, source=2, destination=0, fraction=1.0),
                fwdpy11.copy_individuals(when=5, source=0, destination=2, fraction=1.0),
            ],
            set_deme_sizes=[fwdpy11.SetDemeSize(when=10, deme=2, new_size=10)],
            set_growth_rates=None,
            total_simulation_length=20,
            expected_itree=intervaltree.IntervalTree(
                [
                    intervaltree.Interval(
                        0,
                        2,
                        fwdpy11.discrete_demography._EpochData(0, [0], 10, 1.0),
                    ),
                    intervaltree.Interval(
                        2,
                        21,
                        fwdpy11.discrete_demography._EpochData(0, [0, 1, 2], 30, 1.0),
                    ),
                    intervaltree.Interval(
                        0, 21, fwdpy11.discrete_demography._EpochData(1, [1], 10, 1.0)
                    ),
                    intervaltree.Interval(
                        0, 6, fwdpy11.discrete_demography._EpochData(2, [2], 10, 1.0)
                    ),
                    intervaltree.Interval(
                        6,
                        11,
                        fwdpy11.discrete_demography._EpochData(2, [0, 2], 40, 1.0),
                    ),
                    intervaltree.Interval(
                        11, 21, fwdpy11.discrete_demography._EpochData(2, [2], 10, 1.0)
                    ),
                ]
            ),
        ),
    ],
)
def test_deme_size_history_three_deme_mass_migrations(payload: Payload):
    __detailed_size_history_test_runner(payload)


@pytest.mark.parametrize(
    "payload",
    [
        # One deme is replaced by another at when=1
        # NOTE: this test helped find a bug in our
        # class, but it cannot be run through the back-end
        # b/c there's no ancesty for deme 1
        Payload(
            initial_sizes={0: 10},
            mass_migrations=None,
            set_deme_sizes=[
                fwdpy11.SetDemeSize(deme=0, when=1, new_size=0),
                fwdpy11.SetDemeSize(deme=1, when=1, new_size=10),
            ],
            set_growth_rates=None,
            total_simulation_length=None,
            expected_itree=intervaltree.IntervalTree(
                [
                    intervaltree.Interval(
                        0,
                        2,
                        fwdpy11.discrete_demography._EpochData(0, [0], 10, 1.0),
                    ),
                    intervaltree.Interval(
                        2,
                        np.iinfo(np.uint32).max,
                        fwdpy11.discrete_demography._EpochData(1, None, 10, 1.0),
                    ),
                ]
            ),
        ),
        # NOTE: this is a valid interval tree,
        # but the model cannot be run below b/c there are no migrations
        # to give deme 1 ancestors. Will need to make sure we cover this
        # case when we refactor the DemographyDebugger.
        Payload(
            initial_sizes={0: 10},
            mass_migrations=None,
            set_deme_sizes=[fwdpy11.SetDemeSize(when=0, deme=1, new_size=10)],
            set_growth_rates=None,
            total_simulation_length=20,
            expected_itree=intervaltree.IntervalTree(
                [
                    intervaltree.Interval(
                        0, 21, fwdpy11.discrete_demography._EpochData(0, [0], 10, 1.0)
                    ),
                    intervaltree.Interval(
                        1, 21, fwdpy11.discrete_demography._EpochData(1, None, 10, 1.0)
                    ),
                ]
            ),
        ),
    ],
)
def test_deme_size_history_models_with_no_valid_ancestry(payload: Payload):
    """
    Tests of models that give valid size histories but cannot be run
    though the C++ back-end because some demes do not have valid ancestry.

    These are examples of models that the DemographyDebugger should
    flag as invalid.
    """
    import fwdpy11.discrete_demography

    deme_properties = fwdpy11.discrete_demography._DemeSizeHistory.from_lowlevel(
        initial_sizes=payload.initial_sizes,
        mass_migrations=payload.mass_migrations,
        set_deme_sizes=payload.set_deme_sizes,
        set_growth_rates=payload.set_growth_rates,
        total_simulation_length=payload.total_simulation_length,
    )

    assert deme_properties.epochs == payload.expected_itree


@pytest.mark.parametrize(
    "payload",
    [
        Payload(
            initial_sizes={0: 10},
            mass_migrations=None,
            set_deme_sizes=None,
            # The deme does not have a previous history.
            set_growth_rates=[fwdpy11.SetExponentialGrowth(when=0, deme=1, G=1.1)],
            total_simulation_length=20,
            expected_itree=intervaltree.IntervalTree(),
        ),
        # NOTE: we don't have a test of this
        # case in the existing test suite.
        Payload(
            initial_sizes={0: 10},
            mass_migrations=[
                fwdpy11.move_individuals(
                    when=10, source=0, destination=1, fraction=1.0
                ),
                fwdpy11.move_individuals(
                    when=10, source=1, destination=2, fraction=1.0
                ),
            ],
            set_deme_sizes=None,
            set_growth_rates=None,
            total_simulation_length=20,
            expected_itree=intervaltree.IntervalTree(),
        ),
        Payload(
            initial_sizes={0: 10},
            mass_migrations=[
                fwdpy11.copy_individuals(
                    when=10, source=0, destination=1, fraction=1.0
                ),
                fwdpy11.copy_individuals(
                    when=10, source=1, destination=2, fraction=1.0
                ),
            ],
            set_deme_sizes=None,
            set_growth_rates=None,
            total_simulation_length=20,
            expected_itree=intervaltree.IntervalTree(),
        ),
    ],
)
def test_deme_size_history_bad_models(payload: Payload):
    import fwdpy11.discrete_demography

    with pytest.raises(ValueError):
        _ = fwdpy11.discrete_demography._DemeSizeHistory.from_lowlevel(
            initial_sizes=payload.initial_sizes,
            mass_migrations=payload.mass_migrations,
            set_deme_sizes=payload.set_deme_sizes,
            set_growth_rates=payload.set_growth_rates,
            total_simulation_length=payload.total_simulation_length,
        )


# NOTE: this tests internal implementation details.
# This test should be deleted once a public API
# is in place
def test_invalid_demes_spec_modules():
    for bad in glob.glob("tests/demes-spec/test-cases/invalid/*.yaml"):
        with open(bad, "r") as f:
            yaml = "".join(f.readlines())
            with pytest.raises(fwdpy11.DemographyError):
                _ = fwdpy11._fwdpy11._ForwardDemesGraph(yaml=yaml, burnin=0)


# NOTE: this tests internal implementation details.
# This test should be deleted once a public API
# is in place
def test_valid_demes_spec_modules():
    for good in glob.glob("tests/demes-spec/test-cases/valid/*.yaml"):
        with open(good, "r") as f:
            yaml = "".join(f.readlines())
            _ = fwdpy11._fwdpy11._ForwardDemesGraph(yaml=yaml, burnin=0)


if __name__ == "__main__":
    unittest.main()
