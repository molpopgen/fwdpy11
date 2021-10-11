import copy
import typing
import unittest
from dataclasses import dataclass

import demes
import fwdpy11
import numpy as np
import pytest


def check_valid_demography(cls):
    def _valid_fwdpy11_demography(self):
        try:
            _ = fwdpy11.DemographyDebugger(
                [100] * len(self.demog.metadata["initial_sizes"]), self.demog
            )
        except:
            self.fail("unexpected exception")

    cls.test_validity = _valid_fwdpy11_demography
    return cls


@check_valid_demography
class TestNoEvents(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(description="test demography", time_units="generations")
        self.b.add_deme(name="deme", epochs=[dict(start_size=1000, end_time=0)])

        self.g = self.b.resolve()
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 10)


class TestBadBurnin(unittest.TestCase):
    def test_burnin_inputs(self):
        self.b = demes.Builder(description="test demography", time_units="generations")
        self.b.add_deme(name="deme", epochs=[dict(start_size=1000, end_time=0)])

        self.g = self.b.resolve()
        with pytest.raises(ValueError):
            self.demog = fwdpy11.discrete_demography.from_demes(self.g, 10.0)
        with pytest.raises(ValueError):
            self.demog = fwdpy11.discrete_demography.from_demes(self.g, -1)


@check_valid_demography
class TestLoadGraph(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.g = demes.load("tests/test_demog.yaml")
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 10)


@check_valid_demography
class TestLoadYAML(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.demog = fwdpy11.discrete_demography.from_demes("tests/test_demog.yaml", 10)


@check_valid_demography
class TestTwoEpoch(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(description="test demography", time_units="generations")
        self.b.add_deme(
            name="deme",
            epochs=[
                dict(start_size=1000, end_time=100),
                dict(start_size=2000, end_time=0),
            ],
        )
        self.g = self.b.resolve()
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 10)

    def test_size_change_params(self):
        self.assertEqual(len(self.demog.model.set_deme_sizes), 1)
        self.assertEqual(
            self.demog.model.set_deme_sizes[0].when,
            self.demog.metadata["burnin_time"],
        )
        self.assertEqual(self.demog.model.set_deme_sizes[0].new_size, 2000)


@check_valid_demography
class TestNonGenerationUnits(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(
            description="test demography", time_units="years", generation_time=25
        )
        self.b.add_deme(
            name="Pop",
            epochs=[
                dict(start_size=10000, end_time=25000),
                dict(start_size=1000, end_size=20000, end_time=0),
            ],
        )
        self.g = self.b.resolve()
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 10)

    def test_conversion_to_generations(self):
        self.assertEqual(
            self.demog.model.set_deme_sizes[0].when,
            self.demog.metadata["burnin_time"],
        )
        self.assertEqual(
            self.demog.metadata["total_simulation_length"],
            self.demog.metadata["burnin_time"] + 25000 // 25,
        )
        self.assertEqual(self.demog.model.set_deme_sizes[0].new_size, 1000)


@check_valid_demography
class TestSelfingShift(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(description="test demography", time_units="generations")
        self.b.add_deme(
            name="Selfer",
            epochs=[
                dict(start_size=1000, end_time=1000),
                dict(start_size=1000, end_time=0, selfing_rate=0.2),
            ],
        )
        self.g = self.b.resolve()
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 10)

    def test_selfing_parameters(self):
        pass
        # FIXME: this test relies on internal sorting orders of data,
        # which are not required for correctness
        # self.assertEqual(self.demog.model.set_selfing_rates[0].when, 0)
        # self.assertEqual(self.demog.model.set_selfing_rates[0].S, 0.0)
        # self.assertEqual(
        #     self.demog.model.set_selfing_rates[1].when,
        #     self.demog.metadata["burnin_time"] + 1,
        # )
        # self.assertTrue(self.demog.model.set_selfing_rates[1].S == 0.2)


@check_valid_demography
class TestSelfing(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(description="test demography", time_units="generations")
        self.b.add_deme(
            name="Selfer",
            epochs=[dict(start_size=1000, end_time=0)],
            defaults=dict(epoch=dict(selfing_rate=0.5)),
        )
        self.g = self.b.resolve()
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 10)

    def test_single_pop_selfing(self):
        self.assertTrue(self.demog.model.set_selfing_rates[0].when == 0)
        self.assertTrue(self.demog.model.set_selfing_rates[0].S == 0.5)


@check_valid_demography
class TestSplit(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(description="test demography", time_units="generations")
        self.b.add_deme(name="Ancestor", epochs=[dict(start_size=1000, end_time=200)])
        self.b.add_deme(
            "Deme1", epochs=[dict(start_size=100, end_time=0)], ancestors=["Ancestor"]
        )
        self.b.add_deme(
            "Deme2", epochs=[dict(start_size=100, end_time=0)], ancestors=["Ancestor"]
        )
        self.g = self.b.resolve()
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 10)

    def test_size_changes(self):
        self.assertEqual(len(self.demog.model.set_deme_sizes), 3)

        check_debugger_passes(self.demog)

        # NOTE: commented these out as they test internal sort order and not
        # the model validity
        # self.assertTrue(self.demog.model.set_deme_sizes[0].deme == 1)
        # self.assertTrue(self.demog.model.set_deme_sizes[0].new_size == 100)
        # self.assertTrue(
        #     self.demog.model.set_deme_sizes[0].when
        #     == self.demog.metadata["burnin_time"]
        # )
        # self.assertTrue(self.demog.model.set_deme_sizes[1].deme == 2)
        # self.assertTrue(self.demog.model.set_deme_sizes[1].new_size == 100)
        # self.assertTrue(
        #     self.demog.model.set_deme_sizes[1].when
        #     == self.demog.metadata["burnin_time"]
        # )
        # self.assertTrue(self.demog.model.set_deme_sizes[2].deme == 0)
        # self.assertTrue(self.demog.model.set_deme_sizes[2].new_size == 0)
        # self.assertTrue(
        #     self.demog.model.set_deme_sizes[2].when
        #     == self.demog.metadata["burnin_time"]
        # )


@check_valid_demography
class TestSplitMigration(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(description="test demography", time_units="generations")
        self.b.add_deme(name="Ancestor", epochs=[dict(start_size=1000, end_time=200)])
        self.b.add_deme(
            "Deme1", epochs=[dict(start_size=100, end_time=0)], ancestors=["Ancestor"]
        )
        self.b.add_deme(
            "Deme2", epochs=[dict(start_size=100, end_time=0)], ancestors=["Ancestor"]
        )
        self.b.add_migration(source="Deme1", dest="Deme2", rate=0.01)
        self.b.add_migration(source="Deme2", dest="Deme1", rate=0.02)
        self.g = self.b.resolve()
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 10)


@check_valid_demography
class TestSplitSymmetricMigration(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(description="test demography", time_units="generations")
        self.b.add_deme(name="Ancestor", epochs=[dict(start_size=1000, end_time=200)])
        self.b.add_deme(
            "Deme1", epochs=[dict(start_size=100, end_time=0)], ancestors=["Ancestor"]
        )
        self.b.add_deme(
            "Deme2", epochs=[dict(start_size=100, end_time=0)], ancestors=["Ancestor"]
        )
        self.b.add_migration(demes=["Deme1", "Deme2"], rate=0.01)
        self.g = self.b.resolve()
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 10)


@check_valid_demography
class TestSplitThreeWay(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(description="test demography", time_units="generations")
        self.b.add_deme(name="Ancestor", epochs=[dict(start_size=1000, end_time=200)])
        self.b.add_deme(
            "Deme1", epochs=[dict(start_size=100, end_time=0)], ancestors=["Ancestor"]
        )
        self.b.add_deme(
            "Deme2", epochs=[dict(start_size=100, end_time=0)], ancestors=["Ancestor"]
        )
        self.b.add_deme(
            "Deme3", epochs=[dict(start_size=200, end_time=0)], ancestors=["Ancestor"]
        )
        self.g = self.b.resolve()
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 10)


@check_valid_demography
class TestSplitThreeWayMigration(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(description="test demography", time_units="generations")
        self.b.add_deme(name="Ancestor", epochs=[dict(start_size=1000, end_time=200)])
        self.b.add_deme(
            "Deme1", epochs=[dict(start_size=100, end_time=0)], ancestors=["Ancestor"]
        )
        self.b.add_deme(
            "Deme2", epochs=[dict(start_size=100, end_time=0)], ancestors=["Ancestor"]
        )
        self.b.add_deme(
            "Deme3", epochs=[dict(start_size=200, end_time=0)], ancestors=["Ancestor"]
        )
        self.b.add_migration(demes=["Deme1", "Deme2", "Deme3"], rate=0.1)
        self.g = self.b.resolve()
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 10)


@check_valid_demography
class TestBranch(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(description="test branch", time_units="generations")
        self.b.add_deme(name="Ancestor", epochs=[dict(start_size=1000, end_time=0)])
        self.b.add_deme(
            "Deme1",
            epochs=[dict(start_size=100, end_time=0)],
            ancestors=["Ancestor"],
            start_time=100,
        )
        self.g = self.b.resolve()
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 10)


@check_valid_demography
class TestBranchMigration(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(description="test branch", time_units="generations")
        self.b.add_deme(name="Ancestor", epochs=[dict(start_size=1000, end_time=0)])
        self.b.add_deme(
            "Deme1",
            epochs=[dict(start_size=100, end_time=0)],
            ancestors=["Ancestor"],
            start_time=100,
        )
        self.b.add_migration(source="Ancestor", dest="Deme1", rate=0.01)
        self.b.add_migration(source="Deme1", dest="Ancestor", rate=0.01)
        self.g = self.b.resolve()
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 10)


@check_valid_demography
class TestMultipleBranches(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(description="test branch", time_units="generations")
        self.b.add_deme(name="Ancestor", epochs=[dict(start_size=1000, end_time=0)])
        self.b.add_deme(
            "Deme1",
            epochs=[dict(start_size=100, end_time=0)],
            ancestors=["Ancestor"],
            start_time=100,
        )
        self.b.add_deme(
            "Deme2",
            epochs=[dict(start_size=200, end_time=0)],
            ancestors=["Ancestor"],
            start_time=50,
        )
        self.b.add_deme(
            "Deme3",
            epochs=[dict(start_size=300, end_time=0)],
            ancestors=["Deme1"],
            start_time=20,
        )
        self.g = self.b.resolve()
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 10)


@check_valid_demography
class TestSplitsBranches(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(description="test", time_units="generations")
        self.b.add_deme(name="A", epochs=[dict(start_size=1000, end_time=100)])
        self.b.add_deme(
            name="B",
            epochs=[dict(start_size=1000, end_time=0)],
            ancestors=["A"],
            start_time=200,
        )
        self.b.add_deme(
            name="C", epochs=[dict(start_size=1000, end_time=0)], ancestors=["A"]
        )
        self.b.add_deme(
            name="D", epochs=[dict(start_size=1000, end_time=0)], ancestors=["A"]
        )
        self.g = self.b.resolve()
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 10)


@check_valid_demography
class TestIslandModel(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(description="island", time_units="generations")
        self.b.add_deme(name="Island1", epochs=[dict(start_size=100, end_time=0)])
        self.b.add_deme(name="Island2", epochs=[dict(start_size=200, end_time=0)])
        self.b.add_migration(source="Island1", dest="Island2", rate=0.01)
        self.b.add_migration(source="Island2", dest="Island1", rate=0.02)
        self.g = self.b.resolve()
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 10)

    def test_demog_attributes(self):
        self.assertTrue(
            self.demog.metadata["burnin_time"]
            == sum(self.demog.metadata["initial_sizes"].values()) * 10
        )
        self.assertTrue(len(self.demog.model.set_migration_rates) == 0)


@check_valid_demography
class TestIslandModelRateChange(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(description="island", time_units="generations")
        self.b.add_deme(name="Island1", epochs=[dict(start_size=100, end_time=0)])
        self.b.add_deme(name="Island2", epochs=[dict(start_size=200, end_time=0)])
        self.b.add_migration(source="Island1", dest="Island2", rate=0.01, end_time=500)
        self.b.add_migration(source="Island2", dest="Island1", rate=0.02)
        self.b.add_migration(
            source="Island1", dest="Island2", rate=0.05, start_time=500
        )
        self.g = self.b.resolve()
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 10)

    def test_burnin_time(self):
        self.assertTrue(
            self.demog.metadata["burnin_time"]
            == sum(self.demog.metadata["initial_sizes"].values()) * 10
        )

    def test_num_mig_rate_changes(self):
        self.assertTrue(len(self.demog.model.set_migration_rates) == 1)

    def test_total_sim_length(self):
        self.assertEqual(
            self.demog.metadata["total_simulation_length"],
            self.demog.metadata["burnin_time"] + 500,
        )


@check_valid_demography
class TestTwoPopMerger(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(
            description="split then merger", time_units="generations"
        )
        self.b.add_deme(name="Ancestral", epochs=[dict(start_size=1000, end_time=1000)])
        self.b.add_deme(
            name="Parent1",
            epochs=[dict(start_size=500, end_time=500)],
            ancestors=["Ancestral"],
        )
        self.b.add_deme(
            name="Parent2",
            epochs=[dict(start_size=500, end_time=500)],
            ancestors=["Ancestral"],
        )
        self.b.add_deme(
            name="Child",
            epochs=[dict(start_size=1000, end_time=0)],
            ancestors=["Parent1", "Parent2"],
            proportions=[0.5, 0.5],
            start_time=500,
        )
        self.g = self.b.resolve()
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 10)

    def test_total_sim_length(self):
        self.assertEqual(
            self.demog.metadata["total_simulation_length"],
            self.demog.metadata["burnin_time"] + 1000,
        )

    def test_num_size_changes(self):
        self.assertEqual(len(self.demog.model.set_deme_sizes), 6)


@check_valid_demography
class TestFourWayMerger(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(
            description="split then merger", time_units="generations"
        )
        self.b.add_deme(name="Ancestral", epochs=[dict(start_size=1000, end_time=1000)])
        self.b.add_deme(
            name="A",
            epochs=[dict(start_size=500, end_time=700)],
            ancestors=["Ancestral"],
        )
        self.b.add_deme(
            name="B",
            epochs=[dict(start_size=500, end_time=500)],
            ancestors=["Ancestral"],
        )
        self.b.add_deme(
            name="Parent1", epochs=[dict(start_size=200, end_time=100)], ancestors=["A"]
        )
        self.b.add_deme(
            name="Parent2", epochs=[dict(start_size=300, end_time=100)], ancestors=["A"]
        )
        self.b.add_deme(
            name="Parent3", epochs=[dict(start_size=100, end_time=100)], ancestors=["B"]
        )
        self.b.add_deme(
            name="Parent4", epochs=[dict(start_size=400, end_time=100)], ancestors=["B"]
        )
        self.b.add_deme(
            name="Child",
            epochs=[dict(start_size=1000, end_time=0)],
            ancestors=["Parent1", "Parent2", "Parent3", "Parent4"],
            proportions=[0.2, 0.3, 0.1, 0.4],
            start_time=100,
        )
        self.g = self.b.resolve()
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 10)

    def test_total_sim_length(self):
        self.assertEqual(
            self.demog.metadata["total_simulation_length"],
            self.demog.metadata["burnin_time"] + 1000,
        )

    def test_num_size_changes(self):
        self.assertTrue(len(self.demog.model.set_deme_sizes), 14)


@check_valid_demography
class TestPulseMigration(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(description="pulse", time_units="generations")
        self.b.add_deme(name="deme1", epochs=[dict(start_size=100, end_time=0)])
        self.b.add_deme(name="deme2", epochs=[dict(start_size=100, end_time=0)])
        self.b.add_pulse(source="deme1", dest="deme2", time=100, proportion=0.2)
        self.g = self.b.resolve()
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 10)

    def test_total_sim_length(self):
        self.assertEqual(
            self.demog.metadata["total_simulation_length"],
            self.demog.metadata["burnin_time"] + 100,
        )

    def test_pulse_migration_matrix(self):
        self.assertTrue(len(self.demog.model.set_migration_rates) == 2)
        self.assertEqual(
            self.demog.model.set_migration_rates[0].when,
            self.demog.metadata["burnin_time"],
        )
        self.assertEqual(
            self.demog.model.set_migration_rates[1].when,
            self.demog.metadata["burnin_time"] + 1,
        )
        self.assertTrue(self.demog.model.set_migration_rates[0].deme == 1)
        self.assertTrue(self.demog.model.set_migration_rates[1].deme == 1)
        self.assertTrue(
            np.all(self.demog.model.set_migration_rates[0].migrates == [0.2, 0.8])
        )
        self.assertTrue(
            np.all(self.demog.model.set_migration_rates[1].migrates == [0, 1])
        )


def check_debugger_passes(demog):
    try:
        _ = fwdpy11.DemographyDebugger(
            [100] * len(demog.metadata["initial_sizes"]), demog
        )
    except:
        raise "unexpected exception"


def evolve_demes_model(demog, initial_sizes) -> fwdpy11.DiploidPopulation:
    pdict = {
        "rates": (0.0, 0.0, 0.0),
        "gvalue": fwdpy11.Multiplicative(2.0),
        "demography": demog,
        "simlen": demog.metadata["total_simulation_length"],
    }
    params = fwdpy11.ModelParams(**pdict)
    rng = fwdpy11.GSLrng(100)
    pop = fwdpy11.DiploidPopulation(initial_sizes, 1.0)
    fwdpy11.evolvets(rng, pop, params, 100)
    return pop


def three_way_continuous_migration(
    goes_extinct: typing.Optional[typing.List[str]] = None,
):
    b = demes.Builder(description="many migrations", time_units="generations")
    for deme in ["A", "B", "C"]:
        if goes_extinct is not None and deme in goes_extinct:
            b.add_deme(name=deme, epochs=[dict(start_size=100, end_time=10)])
        else:
            b.add_deme(name=deme, epochs=[dict(start_size=100, end_time=0)])
    b.add_migration(demes=["A", "B", "C"], rate=0.1)
    g = b.resolve()
    return fwdpy11.discrete_demography.from_demes(g, 1)


def three_way_continuous_migration_pairwise(
    goes_extinct: typing.Optional[typing.List[str]] = None,
):
    b = demes.Builder(description="many migrations", time_units="generations")
    for deme in ["A", "B", "C"]:
        if goes_extinct is not None and deme in goes_extinct:
            b.add_deme(name=deme, epochs=[dict(start_size=100, end_time=10)])
        else:
            b.add_deme(name=deme, epochs=[dict(start_size=100, end_time=0)])
    b.add_migration(demes=["A", "B"], rate=0.1)
    b.add_migration(demes=["A", "C"], rate=0.1)
    b.add_migration(demes=["B", "C"], rate=0.1)
    g = b.resolve()
    return fwdpy11.discrete_demography.from_demes(g, 1)


@pytest.mark.parametrize(
    "demog",
    [
        three_way_continuous_migration(),
        three_way_continuous_migration(["A"]),
        three_way_continuous_migration(["B"]),
        three_way_continuous_migration(["C"]),
        three_way_continuous_migration(["A", "B"]),
        three_way_continuous_migration(["A", "C"]),
        three_way_continuous_migration(["B", "C"]),
        three_way_continuous_migration_pairwise(),
        three_way_continuous_migration_pairwise(["A"]),
        three_way_continuous_migration_pairwise(["B"]),
        three_way_continuous_migration_pairwise(["C"]),
        three_way_continuous_migration_pairwise(["A", "B"]),
        three_way_continuous_migration_pairwise(["A", "C"]),
        three_way_continuous_migration_pairwise(["B", "C"]),
    ],
)
def test_three_way_continuous_migration_pairwise(demog):
    check_debugger_passes(demog)
    assert np.all(
        demog.model.migmatrix.M
        == np.array([[0.8, 0.1, 0.1], [0.1, 0.8, 0.1], [0.1, 0.1, 0.8]])
    )


@pytest.mark.parametrize(
    "demog",
    [
        three_way_continuous_migration(),
        three_way_continuous_migration(["A"]),
        three_way_continuous_migration(["B"]),
        three_way_continuous_migration(["C"]),
        three_way_continuous_migration(["A", "B"]),
        three_way_continuous_migration(["A", "C"]),
        three_way_continuous_migration(["B", "C"]),
        three_way_continuous_migration_pairwise(),
        three_way_continuous_migration_pairwise(["A"]),
        three_way_continuous_migration_pairwise(["B"]),
        three_way_continuous_migration_pairwise(["C"]),
        three_way_continuous_migration_pairwise(["A", "B"]),
        three_way_continuous_migration_pairwise(["A", "C"]),
        three_way_continuous_migration_pairwise(["B", "C"]),
    ],
)
def test_evolve_three_way_continuous_migration_pairwise(demog):
    pop = evolve_demes_model(demog, [100] * 3)
    assert pop.generation == demog.metadata["total_simulation_length"]


def multiple_migrations_delayed():
    b = demes.Builder(description="many migrations", time_units="generations")
    b.add_deme(name="A", epochs=[dict(start_size=100, end_time=0)])
    b.add_deme(name="B", epochs=[dict(start_size=100, end_time=0)])
    b.add_deme(name="C", epochs=[dict(start_size=100, end_time=0)])
    b.add_migration(demes=["A", "B"], rate=0.1)
    b.add_migration(demes=["B", "C"], rate=0.1)
    b.add_migration(demes=["A", "C"], rate=0.1, start_time=100)
    g = b.resolve()
    return fwdpy11.discrete_demography.from_demes(g, 1)


@pytest.mark.parametrize("demog", [multiple_migrations_delayed()])
def test_multiple_migrations_delayed(demog):
    check_debugger_passes(demog)
    assert np.all(
        demog.model.migmatrix.M
        == np.array([[0.9, 0.1, 0], [0.1, 0.8, 0.1], [0, 0.1, 0.9]])
    )
    assert len(demog.model.set_migration_rates) == 2

    for set_mig in demog.model.set_migration_rates:
        assert set_mig.deme in [0, 2]
        if set_mig.deme == 0:
            assert np.all(set_mig.migrates == np.array([0.8, 0.1, 0.1]))
        elif set_mig.deme == 2:
            assert np.all(set_mig.migrates == np.array([0.1, 0.1, 0.8]))


def splits_with_migrations():
    b = demes.Builder(description="splits with migration", time_units="generations")
    b.add_deme(name="A", epochs=[dict(start_size=100, end_time=100)])
    b.add_deme(name="B", ancestors=["A"], epochs=[dict(start_size=100, end_time=0)])
    b.add_deme(name="C", ancestors=["A"], epochs=[dict(start_size=100, end_time=50)])
    b.add_deme(name="D", ancestors=["C"], epochs=[dict(start_size=100, end_time=0)])
    b.add_deme(name="E", ancestors=["C"], epochs=[dict(start_size=100, end_time=0)])
    b.add_migration(demes=["B", "C"], rate=0.1)
    b.add_migration(demes=["B", "D"], rate=0.1)
    b.add_migration(demes=["B", "E"], rate=0.1)
    b.add_migration(demes=["D", "E"], rate=0.1)
    g = b.resolve()
    return fwdpy11.discrete_demography.from_demes(g, 1)


@pytest.mark.parametrize("demog", [splits_with_migrations()])
def test_splits_with_migrations(demog):
    check_debugger_passes(demog)
    M = np.zeros(25).reshape(5, 5)
    M[0, 0] = 1
    assert np.all(demog.model.migmatrix.M == M)
    num_splits = 2
    len_set_migs = 3 * num_splits + 2 + 3
    assert len(demog.model.set_migration_rates) == len_set_migs


def yaml_migration_1():
    return """
description: test model
time_units: generations
demes:
  - name: A
    epochs:
      - end_time: 200
        start_size: 100
  - name: B
    ancestors: [A]
    epochs:
      - end_time: 100
        start_size: 100
  - name: C
    ancestors: [A]
    epochs:
    - end_time: 20
      start_size: 100
  - name: D
    ancestors: [B]
    epochs:
      - end_time: 2
        start_size: 100
      - end_time: 0
        start_size: 10
  - name: E
    ancestors: [B]
    epochs:
      - end_time: 2
        start_size: 100
      - end_time: 0
        start_size: 30
migrations:
  - demes: [C, D]
    rate: 6.228e-5
  - demes: [D, E]
    rate: 4.14e-5
"""


def yaml_migration_2():
    """
    Same model as previous yaml with deme name D changed to F to switch sort order.
    """
    return """description: test model
time_units: generations
demes:
  - name: A
    epochs:
      - end_time: 200
        start_size: 100
  - name: B
    ancestors: [A]
    epochs:
      - end_time: 100
        start_size: 100
  - name: C
    ancestors: [A]
    epochs:
    - end_time: 20
      start_size: 100
  - name: F
    ancestors: [B]
    epochs:
      - end_time: 2
        start_size: 100
      - end_time: 0
        start_size: 10
  - name: E
    ancestors: [B]
    epochs:
      - end_time: 2
        start_size: 100
      - end_time: 0
        start_size: 30
migrations:
  - demes: [C, F]
    rate: 6.228e-5
  - demes: [F, E]
    rate: 4.14e-5
"""


@pytest.mark.parametrize("data", [yaml_migration_1(), yaml_migration_2()])
def test_yamls_with_migration(data):
    g = demes.loads(data)
    demog = fwdpy11.discrete_demography.from_demes(g, 1)
    check_debugger_passes(demog)


def test_split_model_population_size_history(two_deme_split_with_ancestral_size_change):
    """
    This is a detailed test of the complete size history
    of a model that we run all the way through.
    """
    model = fwdpy11.discrete_demography.from_demes(
        two_deme_split_with_ancestral_size_change, burnin=1
    )

    @dataclass
    class DemeSizeAtTime:
        when: int
        size: int

    class DemeSizes(object):
        def __init__(self):
            self.sizes = dict()

        def __call__(self, pop, _):
            for key, value in pop.deme_sizes(as_dict=True).items():
                if key not in self.sizes:
                    self.sizes[key] = [DemeSizeAtTime(when=pop.generation, size=value)]
                else:
                    self.sizes[key].append(
                        DemeSizeAtTime(when=pop.generation, size=value)
                    )

    pdict = {
        "gvalue": fwdpy11.Multiplicative(2.0),
        "rates": (0, 0, 0),
        "demography": model,
        "simlen": model.metadata["total_simulation_length"],
    }
    params = fwdpy11.ModelParams(**pdict)
    pop = fwdpy11.DiploidPopulation(100, 1.0)
    rng = fwdpy11.GSLrng(90210)
    recorder = DemeSizes()
    fwdpy11.evolvets(rng, pop, params, 100, recorder=recorder)

    # The ancestral deme exists until generation 110,
    # and we only see offspring from birth time 1 on.
    assert [i.when for i in recorder.sizes[0]] == [i for i in range(1, 111)]
    # The daughter demes are seen from 110 till the end
    for deme in [1, 2]:
        assert [i.when for i in recorder.sizes[deme]] == [
            i for i in range(111, model.metadata["total_simulation_length"] + 1)
        ]
    # initial daughter deme sizes
    g1 = 2 ** (1 / 10) - 1
    g2 = 4 ** (1 / 10) - 1
    assert recorder.sizes[1][0].size == int(np.rint(250 * (1 + g1)))
    assert recorder.sizes[2][0].size == int(np.rint(50 * (1 + g2)))
    # final daughter deme sizes
    assert recorder.sizes[1][-1].size == 500
    assert recorder.sizes[2][-1].size == 200

    # At generation 100, the ancestral pop size changed from 100
    # to 200
    for i in recorder.sizes[0]:
        if i.when < 101:
            assert i.size == 100, f"{i}"
        else:
            assert i.size == 200, f"{i}"


@pytest.mark.parametrize("when", [i for i in range(75, 120)])
def test_evolve_population_in_two_stages(
    when, two_deme_split_with_ancestral_size_change
):
    model = fwdpy11.discrete_demography.from_demes(
        two_deme_split_with_ancestral_size_change, burnin=1
    )
    pdict = {
        "gvalue": fwdpy11.Multiplicative(2.0),
        "rates": (0, 0, 0),
        "demography": model,
        "simlen": when,
    }
    params = fwdpy11.ModelParams(**pdict)
    pop = fwdpy11.DiploidPopulation(100, 1.0)
    rng = fwdpy11.GSLrng(90210)
    fwdpy11.evolvets(rng, pop, params, 100)

    pdict["simlen"] = model.metadata["total_simulation_length"] - when
    params = fwdpy11.ModelParams(**pdict)

    fwdpy11.evolvets(rng, pop, params, 100, check_demographic_event_timings=False)

    counts = np.unique(np.array(pop.diploid_metadata)["deme"], return_counts=True)
    assert counts[1][0] == 500
    assert counts[1][1] == 200


@pytest.mark.parametrize("when", [i for i in range(75, 120)])
def test_evolve_population_in_two_stages_with_deepcopy(
    when, two_deme_split_with_ancestral_size_change
):
    model = fwdpy11.discrete_demography.from_demes(
        two_deme_split_with_ancestral_size_change, burnin=1
    )
    pdict = {
        "gvalue": fwdpy11.Multiplicative(2.0),
        "rates": (0, 0, 0),
        "demography": model,
        "simlen": when,
    }
    params = fwdpy11.ModelParams(**pdict)
    pop = fwdpy11.DiploidPopulation(100, 1.0)
    rng = fwdpy11.GSLrng(90210)
    fwdpy11.evolvets(rng, pop, params, 100)

    pdict2 = copy.deepcopy(pdict)
    pdict2["simlen"] = model.metadata["total_simulation_length"] - when
    params2 = fwdpy11.ModelParams(**pdict2)

    fwdpy11.evolvets(rng, pop, params2, 100, check_demographic_event_timings=False)

    counts = np.unique(np.array(pop.diploid_metadata)["deme"], return_counts=True)
    assert counts[1][0] == 500, f"{counts}"
    assert counts[1][1] == 200, f"{counts}"


# NOTE: update this test to have burnin=0 once GittHub issue 776
# is fixed
# NOTE: models like this are important to applications where
# the ancestor is simulated with msprime, and the split part
# happens in fwdpy11
@pytest.mark.parametrize("burnin", [1])
def test_evolve_demes_model_starting_with_two_pops_and_no_ancestry(
    burnin,
    start_demes_model_with_two_pops,
):
    model = fwdpy11.discrete_demography.from_demes(
        start_demes_model_with_two_pops, burnin=burnin
    )
    pdict = {
        "gvalue": fwdpy11.Multiplicative(2.0),
        "rates": (0, 0, 0),
        "demography": model,
        "simlen": model.metadata["total_simulation_length"],
    }
    params = fwdpy11.ModelParams(**pdict)
    initial_sizes = [i for i in model.metadata["initial_sizes"].values()]
    pop = fwdpy11.DiploidPopulation(initial_sizes, 1.0)
    for key, value in pop.deme_sizes(as_dict=True).items():
        if key == 0:
            assert value == 100
        elif key == 1:
            assert value == 75
        else:
            raise RuntimeError("unexpected key")
    rng = fwdpy11.GSLrng(101)
    fwdpy11.evolvets(rng, pop, params, 100)
    for key, value in pop.deme_sizes(as_dict=True).items():
        if key == 0:
            assert value == 200
        elif key == 1:
            assert value == 300
        else:
            raise RuntimeError("unexpected key")


def test_four_deme_split_model():
    """
    Tests primarily for lack of temporal
    overlap b/w ancestral/derived demes.
    See GitHub issue #814
    """
    yaml = """time_units: generations
demes:
- name: A
  epochs:
    - {end_time: 45, start_size: 10}
- name: B
  ancestors: [A] 
  epochs:
    - {end_time: 30, start_size: 20}
- name: C
  ancestors: [B] 
  epochs:
    - {end_time: 15, start_size: 30}
- name: D
  ancestors: [C] 
  epochs:
    - {end_time: 0, start_size: 40}
"""
    burnin = 1
    g = demes.loads(yaml)
    model = fwdpy11.discrete_demography.from_demes(g, burnin=burnin)

    @dataclass
    class DemeSizeAtTime:
        when: int
        size: int

    class DemeSizes(object):
        def __init__(self):
            self.sizes = dict()

        def __call__(self, pop, _):
            deme_sizes = pop.deme_sizes()
            assert len(deme_sizes[0]) == 1
            for key, value in pop.deme_sizes(as_dict=True).items():
                if key not in self.sizes:
                    self.sizes[key] = [DemeSizeAtTime(when=pop.generation, size=value)]
                else:
                    self.sizes[key].append(
                        DemeSizeAtTime(when=pop.generation, size=value)
                    )

    pdict = {
        "gvalue": fwdpy11.Multiplicative(2.0),
        "rates": (0, 0, 0),
        "demography": model,
        "simlen": model.metadata["total_simulation_length"],
    }
    params = fwdpy11.ModelParams(**pdict)
    initial_size = g["A"].epochs[0].start_size
    pop = fwdpy11.DiploidPopulation(initial_size, 1.0)
    rng = fwdpy11.GSLrng(90210)
    recorder = DemeSizes()
    fwdpy11.evolvets(rng, pop, params, 100, recorder=recorder)

    for deme in recorder.sizes:
        if deme > 0:
            assert len(recorder.sizes[deme]) == 15
        else:
            assert len(recorder.sizes[deme]) == burnin * initial_size


def no_demography_no_burnin():
    yaml = """time_units: generations
demes:
- name: A
  epochs:
  - {end_time: 0, start_size: 10}
"""
    burnin = 0
    g = demes.loads(yaml)
    model = fwdpy11.discrete_demography.from_demes(g, burnin=burnin)

    @dataclass
    class DemeSizeAtTime:
        when: int
        size: int

    class DemeSizes(object):
        def __init__(self):
            self.sizes = dict()

        def __call__(self, pop, _):
            deme_sizes = pop.deme_sizes()
            assert len(deme_sizes[0]) == 1
            for key, value in pop.deme_sizes(as_dict=True).items():
                if key not in self.sizes:
                    self.sizes[key] = [DemeSizeAtTime(when=pop.generation, size=value)]
                else:
                    self.sizes[key].append(
                        DemeSizeAtTime(when=pop.generation, size=value)
                    )

    pdict = {
        "gvalue": fwdpy11.Multiplicative(2.0),
        "rates": (0, 0, 0),
        "demography": model,
        "simlen": model.metadata["total_simulation_length"],
    }
    params = fwdpy11.ModelParams(**pdict)
    initial_size = g["A"].epochs[0].start_size
    pop = fwdpy11.DiploidPopulation(initial_size, 1.0)
    rng = fwdpy11.GSLrng(90210)
    recorder = DemeSizes()
    fwdpy11.evolvets(rng, pop, params, 100, recorder=recorder)

    assert len(recorder.sizes[0]) == 1
    assert recorder.sizes[0][0].when == 1
    assert recorder.sizes[0][0].size == 10


# The following fixtures are demes models
# with an event that should affect the first generation
# of the forward sim


@dataclass
class PedigreeRow:
    when: int
    deme: int
    label: int
    parents: list
    parent_deme: int

    def is_migrant(self) -> bool:
        return self.deme != self.parent_deme

    def is_selfed(self) -> bool:
        return self.parents[0] == self.parents[1]


@dataclass
class CrudePedigree:
    records: typing.List[PedigreeRow]
    individual_to_deme: list
    individual_record: list

    def __call__(self, pop):
        temp, temp2 = [], []
        for i, md in enumerate(pop.diploid_metadata):
            assert i == md.label

            row = PedigreeRow(pop.generation, md.deme, i, [-1, -1], md.deme)

            if len(self.individual_to_deme) > 0:
                row.parents = [self.individual_record[i] for i in md.parents]
                assert (
                    self.individual_to_deme[md.parents[0]]
                    == self.individual_to_deme[md.parents[1]]
                )
                row.parent_deme = self.individual_to_deme[md.parents[0]]
            else:
                assert pop.generation == 0

            self.records.append(row)
            temp.append(md.deme)
            temp2.append(len(self.records) - 1)

        self.individual_to_deme = temp
        self.individual_record = temp2


@dataclass
class TrackedInfo:
    when: int
    deme: int
    nselfed: int
    demesize: int
    popsize: int


@dataclass
class InfoTracker:
    data: typing.List[TrackedInfo]
    pedigree: CrudePedigree

    def __call__(self, pop, _):
        deme_sizes = pop.deme_sizes(as_dict=True)
        selfed_by_deme = {}
        for i in deme_sizes.keys():
            selfed_by_deme[i] = 0

        for md in pop.diploid_metadata:
            if md.parents[0] == md.parents[1]:
                selfed_by_deme[md.deme] += 1

        for deme, size in deme_sizes.items():
            self.data.append(
                TrackedInfo(pop.generation, deme, selfed_by_deme[deme], size, pop.N)
            )

        self.pedigree(pop)


def build_tracker() -> InfoTracker:
    return InfoTracker([], CrudePedigree([], [], []))


def base_demes_model():
    yaml = """
    time_units: generations
    demes:
    - name: ancestor
      epochs:
        - {end_time: 2, start_size: 100}"""
    return yaml


def set_selfing_rate_generation_1():
    yaml = """
    - name: derived
      ancestors: [ancestor]
      epochs:
        - {end_time: 0, selfing_rate: 1.0, start_size: 100}
    """

    def validate(tracker: InfoTracker):
        for i in tracker.data:
            if i.deme == 0:
                assert i.when == 0
            else:
                assert i.when > 0
                if i.when == 2:
                    assert i.nselfed == i.demesize
                else:
                    # Not rigorous, but the probability
                    # of this being true is quite small
                    assert i.nselfed < i.demesize

        for i in tracker.pedigree.records:
            if i.deme == 1 and i.when == 1:
                assert i.parent_deme == 0

    return base_demes_model() + yaml, build_tracker(), validate


def set_deme_size_generation_1():
    yaml = """
    - name: derived
      ancestors: [ancestor]
      epochs:
        - {end_time: 0, start_size: 66}
    """

    def validate(tracker: InfoTracker):
        for i in tracker.data:
            if i.deme < 1:
                assert i.demesize == 100
                assert i.when == 0
            else:
                assert i.demesize == 66
                assert i.when > 0

    return base_demes_model() + yaml, build_tracker(), validate


def set_growth_rate_generation_1():
    yaml = """
    - name: derived
      ancestors: [ancestor]
      epochs:
        - {end_time: 0, start_size: 50, end_size: 250}
    """

    def validate(tracker: InfoTracker):
        assert tracker.data[-1].demesize == 250
        # Check the initial size of deme 1
        g = 5 ** (1 / 2) - 1
        for i in tracker.data:
            if i.deme == 1:
                assert i.when > 0
                assert i.demesize == np.rint(50 * (1 + g))
                break
            else:
                assert i.when == 0

    return base_demes_model() + yaml, build_tracker(), validate


def simple_split_generation_1():
    yaml = """
    - name: derived0
      ancestors: [ancestor]
      epochs:
        - {end_time: 0, start_size: 50}
    - name: derived1
      ancestors: [ancestor]
      epochs:
        - {end_time: 0, start_size: 50}
    """

    def validate(tracker: InfoTracker):
        for i in tracker.data:
            if i.deme < 1:
                assert i.when < 1
                assert i.demesize == 100
            else:
                assert i.when > 0
                assert i.demesize == 50

    return base_demes_model() + yaml, build_tracker(), validate


def simple_split_with_migration_generation_1():
    yaml = """
    - name: derived0
      ancestors: [ancestor]
      epochs:
        - {end_time: 0, start_size: 50}
    - name: derived1
      ancestors: [ancestor]
      epochs:
        - {end_time: 0, start_size: 50}
    migrations:
      - source: derived0
        dest: derived1
        rate: 1.0
    """

    def validate(tracker: InfoTracker):
        for i in tracker.data:
            if i.deme < 1:
                assert i.when < 1
                assert i.demesize == 100
            else:
                assert i.when > 0
                assert i.demesize == 50

        for i in tracker.pedigree.records:
            if i.deme == 2:
                assert i.is_migrant()
                if i.when > 1:
                    assert i.parent_deme == 1
                else:
                    assert i.parent_deme == 0
            elif i.deme == 1:
                if i.when == 2:
                    assert not i.is_migrant()
                else:
                    assert i.is_migrant()

        # All parents of demes 1 and 2 must be from 0 at generation 1
        for i in tracker.pedigree.records:
            if i.when == 1:
                assert i.deme == 1 or i.deme == 2
                assert i.parent_deme == 0

        deme1parents, deme2parents = [], []
        for i in tracker.pedigree.records:
            if i.when == 2:
                if i.deme == 1:
                    deme1parents.append(i.parent_deme)
                elif i.deme == 2:
                    deme2parents.append(i.parent_deme)

        assert all([i == 1 for i in deme1parents])
        assert not any([i == 2 for i in deme2parents])

    return base_demes_model() + yaml, build_tracker(), validate


def simple_split_with_growth_generation_1():
    yaml = """
    - name: derived0
      ancestors: [ancestor]
      epochs:
        - {end_time: 0, start_size: 50, end_size: 100}
    - name: derived1
      ancestors: [ancestor]
      epochs:
        - {end_time: 0, start_size: 75, end_size: 1000}
    """

    def validate(tracker: InfoTracker):
        for i in tracker.data:
            if i.deme < 1:
                assert i.when == 0
            else:
                assert i.when > 0

        for i in tracker.data[::-1]:
            if i.deme == 1:
                assert i.demesize == 100
                break

        for i in tracker.data[::-1]:
            if i.deme == 2:
                assert i.demesize == 1000
                break

        g = (100 / 50) ** (1 / 2) - 1
        for i in tracker.data:
            if i.deme == 1:
                assert i.demesize == np.rint(50 * (1 + g))
                break

        g = (1000 / 75) ** (1 / 2) - 1
        for i in tracker.data:
            if i.deme == 2:
                assert i.demesize == np.rint(75 * (1 + g))
                break

        # All parents of demes 1 and 2 must be from 0 at generation 1
        for i in tracker.pedigree.records:
            if i.when == 1:
                assert i.deme == 1 or i.deme == 2
                assert i.parent_deme == 0

        deme1parents, deme2parents = [], []
        for i in tracker.pedigree.records:
            if i.when == 2:
                if i.deme == 1:
                    deme1parents.append(i.parent_deme)
                elif i.deme == 2:
                    deme2parents.append(i.parent_deme)

        assert all([i == 1 for i in deme1parents])
        assert all([i == 2 for i in deme2parents])

    return base_demes_model() + yaml, build_tracker(), validate


def simple_split_with_growth_and_migration_generation_1():
    yaml = """
    - name: derived0
      ancestors: [ancestor]
      epochs:
        - {end_time: 0, start_size: 50, end_size: 100}
    - name: derived1
      ancestors: [ancestor]
      epochs:
        - {end_time: 0, start_size: 75, end_size: 1000}
    migrations:
      - source: derived0
        dest: derived1
        rate: 0.01
    """

    def validate(tracker: InfoTracker):
        for i in tracker.data:
            if i.deme < 1:
                assert i.when == 0
            else:
                assert i.when > 0

        for i in tracker.data[::-1]:
            if i.deme == 1:
                assert i.demesize == 100
                break

        for i in tracker.data[::-1]:
            if i.deme == 2:
                assert i.when == 2
                assert i.demesize == 1000
                break

        g = (100 / 50) ** (1 / 2) - 1
        for i in tracker.data:
            if i.deme == 1:
                assert i.when == 1
                assert i.demesize == np.rint(50 * (1 + g))
                break

        g = (1000 / 75) ** (1 / 2) - 1
        for i in tracker.data:
            if i.deme == 2:
                assert i.demesize == np.rint(75 * (1 + g))
                break

        # All parents of demes 1 and 2 must be from 0 at generation 1
        for i in tracker.pedigree.records:
            if i.when == 1:
                assert i.deme == 1 or i.deme == 2
                assert i.parent_deme == 0

        deme1parents, deme2parents = [], []
        for i in tracker.pedigree.records:
            if i.when == 2:
                if i.deme == 1:
                    deme1parents.append(i.parent_deme)
                elif i.deme == 2:
                    deme2parents.append(i.parent_deme)

        assert all([i == 1 for i in deme1parents])
        assert any([i == 2 for i in deme2parents])

    return base_demes_model() + yaml, build_tracker(), validate


def one_generation_hop_generation_1():
    yaml = """
    - name: transient
      ancestors: [ancestor]
      epochs:
        - {end_time: 1, start_size: 100}
    - name: final
      ancestors: [transient]
      epochs:
        - {end_time: 0, start_size: 200}
    """

    def validate(tracker: InfoTracker):
        demes_seen = {}
        for i in tracker.data:
            demes_seen[i.deme] = True
            if i.deme == 0:
                assert i.when == 0
            elif i.deme == 1:
                assert i.when == 1
                assert i.demesize == 100
            else:
                assert i.when == 2
                assert i.demesize == 200

        assert 0 in demes_seen
        assert 1 in demes_seen
        assert 2 in demes_seen

    return base_demes_model() + yaml, build_tracker(), validate


@pytest.mark.parametrize(
    "testdata",
    [
        set_selfing_rate_generation_1,
        set_deme_size_generation_1,
        set_growth_rate_generation_1,
        simple_split_generation_1,
        simple_split_with_migration_generation_1,
        simple_split_with_growth_generation_1,
        simple_split_with_growth_and_migration_generation_1,
        one_generation_hop_generation_1,
    ],
)
@pytest.mark.parametrize("burnin", [0, 1])
def test_events_in_generation_one_following_demes_import(testdata, burnin):
    (model, recorder, validate) = testdata()
    g = demes.loads(model)
    demog = fwdpy11.discrete_demography.from_demes(g, burnin=0)
    pdict = {
        "gvalue": fwdpy11.Multiplicative(2.0),
        "rates": (0, 0, 0),
        "demography": demog,
        "simlen": demog.metadata["total_simulation_length"],
    }
    params = fwdpy11.ModelParams(**pdict)

    initial_sizes = []
    for d in sorted(demog.metadata["initial_sizes"].keys()):
        initial_sizes.append(demog.metadata["initial_sizes"][d])
    pop = fwdpy11.DiploidPopulation(initial_sizes, 1.0)
    recorder(pop, None)
    rng = fwdpy11.GSLrng(90210)
    fwdpy11.evolvets(rng, pop, params, 100, recorder=recorder)

    validate(recorder)


@pytest.mark.parametrize(
    "testdata",
    [
        set_selfing_rate_generation_1,
        set_deme_size_generation_1,
        set_growth_rate_generation_1,
        simple_split_generation_1,
        simple_split_with_migration_generation_1,
        simple_split_with_growth_generation_1,
        simple_split_with_growth_and_migration_generation_1,
        one_generation_hop_generation_1,
    ],
)
@pytest.mark.parametrize("burnin", [0, 1])
def test_events_in_generation_one_following_demes_import_start_stop(testdata, burnin):
    (model, recorder, validate) = testdata()
    g = demes.loads(model)
    demog = fwdpy11.discrete_demography.from_demes(g, burnin=0)
    pdict = {
        "gvalue": fwdpy11.Multiplicative(2.0),
        "rates": (0, 0, 0),
        "demography": demog,
        "simlen": 1,
    }
    params = fwdpy11.ModelParams(**pdict)

    initial_sizes = []
    for d in sorted(demog.metadata["initial_sizes"].keys()):
        initial_sizes.append(demog.metadata["initial_sizes"][d])
    pop = fwdpy11.DiploidPopulation(initial_sizes, 1.0)
    recorder(pop, None)
    rng = fwdpy11.GSLrng(90210)
    fwdpy11.evolvets(rng, pop, params, 100, recorder=recorder)

    pdict["simlen"] = demog.metadata["total_simulation_length"] - 1
    fwdpy11.evolvets(
        rng, pop, params, 100, recorder=recorder, check_demographic_event_timings=False
    )

    validate(recorder)
