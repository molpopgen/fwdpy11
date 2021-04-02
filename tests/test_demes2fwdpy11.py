import demes
import numpy as np
import unittest
import pytest
import fwdpy11


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
        self.assertTrue(len(self.demog.model.set_deme_sizes) == 1)
        self.assertTrue(
            self.demog.model.set_deme_sizes[0].when
            == self.demog.metadata["burnin_time"] - 1
        )
        self.assertTrue(self.demog.model.set_deme_sizes[0].new_size == 2000)


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
        self.assertTrue(
            self.demog.model.set_deme_sizes[0].when
            == self.demog.metadata["burnin_time"] - 1
        )
        self.assertTrue(
            self.demog.metadata["total_simulation_length"]
            == self.demog.metadata["burnin_time"] + 25000 // 25
        )
        self.assertTrue(self.demog.model.set_deme_sizes[0].new_size == 1000)


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
        self.assertTrue(self.demog.model.set_selfing_rates[0].when == 0)
        self.assertTrue(self.demog.model.set_selfing_rates[0].S == 0.0)
        self.assertTrue(
            self.demog.model.set_selfing_rates[1].when
            == self.demog.metadata["burnin_time"]
        )
        self.assertTrue(self.demog.model.set_selfing_rates[1].S == 0.2)


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
        self.assertTrue(len(self.demog.model.set_deme_sizes) == 3)
        self.assertTrue(self.demog.model.set_deme_sizes[0].deme == 1)
        self.assertTrue(self.demog.model.set_deme_sizes[0].new_size == 100)
        self.assertTrue(
            self.demog.model.set_deme_sizes[0].when
            == self.demog.metadata["burnin_time"] - 1
        )
        self.assertTrue(self.demog.model.set_deme_sizes[1].deme == 2)
        self.assertTrue(self.demog.model.set_deme_sizes[1].new_size == 100)
        self.assertTrue(
            self.demog.model.set_deme_sizes[1].when
            == self.demog.metadata["burnin_time"] - 1
        )
        self.assertTrue(self.demog.model.set_deme_sizes[2].deme == 0)
        self.assertTrue(self.demog.model.set_deme_sizes[2].new_size == 0)
        self.assertTrue(
            self.demog.model.set_deme_sizes[2].when
            == self.demog.metadata["burnin_time"]
        )


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
        self.assertTrue(
            self.demog.metadata["total_simulation_length"]
            == self.demog.metadata["burnin_time"] + 500
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
        self.assertTrue(
            self.demog.metadata["total_simulation_length"]
            == self.demog.metadata["burnin_time"] + 1000
        )

    def test_num_size_changes(self):
        self.assertTrue(len(self.demog.model.set_deme_sizes) == 6)


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
        self.assertTrue(
            self.demog.metadata["total_simulation_length"]
            == self.demog.metadata["burnin_time"] + 1000
        )

    def test_num_size_changes(self):
        self.assertTrue(len(self.demog.model.set_deme_sizes) == 14)


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
        self.assertTrue(
            self.demog.metadata["total_simulation_length"]
            == self.demog.metadata["burnin_time"] + 100
        )

    def test_pulse_migration_matrix(self):
        self.assertTrue(len(self.demog.model.set_migration_rates) == 2)
        self.assertTrue(
            self.demog.model.set_migration_rates[0].when
            == self.demog.metadata["burnin_time"] - 1
        )
        self.assertTrue(
            self.demog.model.set_migration_rates[1].when
            == self.demog.metadata["burnin_time"]
        )
        self.assertTrue(self.demog.model.set_migration_rates[0].deme == 1)
        self.assertTrue(self.demog.model.set_migration_rates[1].deme == 1)
        self.assertTrue(
            np.all(self.demog.model.set_migration_rates[0].migrates == [0.2, 0.8])
        )
        self.assertTrue(
            np.all(self.demog.model.set_migration_rates[1].migrates == [0, 1])
        )
