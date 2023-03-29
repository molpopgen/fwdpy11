import copy
import typing
import unittest
from dataclasses import dataclass

import demes
import fwdpy11
import numpy as np
import pytest


def run_model_round_trip(cls):
    def _test_evolvets_roundtrip(self):
        _ = evolve_demes_model(self.demog)

    cls.test_evolvets_roundtrip = _test_evolvets_roundtrip
    return cls


def resolve_and_run(b: demes.Builder, burnin: int = 1):
    g = b.resolve()
    demog = fwdpy11.discrete_demography.from_demes(g, burnin)
    evolve_demes_model(demog)


@run_model_round_trip
class TestNoEvents(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(
            description="test demography", time_units="generations")
        self.b.add_deme(name="deme", epochs=[
                        dict(start_size=1000, end_time=0)])

        self.g = self.b.resolve()
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 1)


class TestBadBurnin(unittest.TestCase):
    def test_burnin_inputs(self):
        self.b = demes.Builder(
            description="test demography", time_units="generations")
        self.b.add_deme(name="deme", epochs=[
                        dict(start_size=1000, end_time=0)])

        self.g = self.b.resolve()
        with pytest.raises(ValueError):
            self.demog = fwdpy11.discrete_demography.from_demes(self.g, 10.0)
        with pytest.raises(ValueError):
            self.demog = fwdpy11.discrete_demography.from_demes(self.g, -1)


@run_model_round_trip
class TestLoadGraph(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.g = demes.load("tests/test_demog.yaml")
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 1)


@run_model_round_trip
class TestLoadYAML(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.demog = fwdpy11.discrete_demography.from_demes(
            "tests/test_demog.yaml", 1)


@run_model_round_trip
class TestTwoEpoch(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(
            description="test demography", time_units="generations")
        self.b.add_deme(
            name="deme",
            epochs=[
                dict(start_size=1000, end_time=100),
                dict(start_size=2000, end_time=0),
            ],
        )
        self.g = self.b.resolve()
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 1)


@run_model_round_trip
class TestNonGenerationUnits(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(
            description="test demography", time_units="years", generation_time=25
        )
        self.b.add_deme(
            name="Pop",
            epochs=[
                dict(start_size=1000, end_time=250),
                dict(start_size=100, end_size=200, end_time=0),
            ],
        )
        self.g = self.b.resolve()
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 1)


@run_model_round_trip
class TestSelfingShift(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(
            description="test demography", time_units="generations")
        self.b.add_deme(
            name="Selfer",
            epochs=[
                dict(start_size=1000, end_time=1000),
                dict(start_size=1000, end_time=0, selfing_rate=0.2),
            ],
        )
        self.g = self.b.resolve()
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 1)


@run_model_round_trip
class TestSelfing(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(
            description="test demography", time_units="generations")
        self.b.add_deme(
            name="Selfer",
            epochs=[dict(start_size=1000, end_time=0)],
            defaults=dict(epoch=dict(selfing_rate=0.5)),
        )
        self.g = self.b.resolve()
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 1)


@run_model_round_trip
class TestSplit(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(
            description="test demography", time_units="generations")
        self.b.add_deme(name="Ancestor", epochs=[
                        dict(start_size=1000, end_time=200)])
        self.b.add_deme(
            "Deme1", epochs=[dict(start_size=100, end_time=0)], ancestors=["Ancestor"]
        )
        self.b.add_deme(
            "Deme2", epochs=[dict(start_size=100, end_time=0)], ancestors=["Ancestor"]
        )
        self.g = self.b.resolve()
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 1)


@run_model_round_trip
class TestSplitMigration(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(
            description="test demography", time_units="generations")
        self.b.add_deme(name="Ancestor", epochs=[
                        dict(start_size=1000, end_time=200)])
        self.b.add_deme(
            "Deme1", epochs=[dict(start_size=100, end_time=0)], ancestors=["Ancestor"]
        )
        self.b.add_deme(
            "Deme2", epochs=[dict(start_size=100, end_time=0)], ancestors=["Ancestor"]
        )
        self.b.add_migration(source="Deme1", dest="Deme2", rate=0.01)
        self.b.add_migration(source="Deme2", dest="Deme1", rate=0.02)
        self.g = self.b.resolve()
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 1)


@run_model_round_trip
class TestSplitSymmetricMigration(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(
            description="test demography", time_units="generations")
        self.b.add_deme(name="Ancestor", epochs=[
                        dict(start_size=1000, end_time=200)])
        self.b.add_deme(
            "Deme1", epochs=[dict(start_size=100, end_time=0)], ancestors=["Ancestor"]
        )
        self.b.add_deme(
            "Deme2", epochs=[dict(start_size=100, end_time=0)], ancestors=["Ancestor"]
        )
        self.b.add_migration(demes=["Deme1", "Deme2"], rate=0.01)
        self.g = self.b.resolve()
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 1)


@run_model_round_trip
class TestSplitAsymmetricMigration(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(
            description="test demography", time_units="generations")
        self.b.add_deme(name="Ancestor", epochs=[
                        dict(start_size=1000, end_time=200)])
        self.b.add_deme(
            "Deme1", epochs=[dict(start_size=100, end_time=0)], ancestors=["Ancestor"]
        )
        self.b.add_deme(
            "Deme2", epochs=[dict(start_size=100, end_time=0)], ancestors=["Ancestor"]
        )
        self.b.add_migration(source="Deme1", dest="Deme2", rate=0.01)
        self.b.add_migration(source="Deme2", dest="Deme1", rate=1e-3)
        self.g = self.b.resolve()
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 1)


@run_model_round_trip
class TestSplitThreeWay(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(
            description="test demography", time_units="generations")
        self.b.add_deme(name="Ancestor", epochs=[
                        dict(start_size=1000, end_time=200)])
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
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 1)


@run_model_round_trip
class TestSplitThreeWayMigration(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(
            description="test demography", time_units="generations")
        self.b.add_deme(name="Ancestor", epochs=[
                        dict(start_size=1000, end_time=200)])
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
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 1)


@run_model_round_trip
class TestBranch(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(description="test branch",
                               time_units="generations")
        self.b.add_deme(name="Ancestor", epochs=[
                        dict(start_size=1000, end_time=0)])
        self.b.add_deme(
            "Deme1",
            epochs=[dict(start_size=100, end_time=0)],
            ancestors=["Ancestor"],
            start_time=100,
        )
        self.g = self.b.resolve()
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 1)


@run_model_round_trip
class TestBranchMigration(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(description="test branch",
                               time_units="generations")
        self.b.add_deme(name="Ancestor", epochs=[
                        dict(start_size=1000, end_time=0)])
        self.b.add_deme(
            "Deme1",
            epochs=[dict(start_size=100, end_time=0)],
            ancestors=["Ancestor"],
            start_time=100,
        )
        self.b.add_migration(source="Ancestor", dest="Deme1", rate=0.01)
        self.b.add_migration(source="Deme1", dest="Ancestor", rate=0.005)
        self.g = self.b.resolve()
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 1)


@pytest.mark.parametrize("second_ancestor", ["Ancestor", "Deme1"])
@pytest.mark.parametrize("first_start_time", [100])
@pytest.mark.parametrize("second_start_time", [99])
@pytest.mark.parametrize(
    "migrations",
    [
        [],
        [
            {"source": "Ancestor", "dest": "Deme1", "rate": 0.01},
        ],
        [
            {"source": "Ancestor", "dest": "Deme2", "rate": 0.01},
        ],
        [
            {"source": "Ancestor", "dest": "Deme1", "rate": 0.01},
            {"source": "Deme1", "dest": "Ancestor", "rate": 0.005},
        ],
        [
            {"source": "Deme2", "dest": "Deme1", "rate": 0.01},
            {"source": "Deme1", "dest": "Deme2", "rate": 0.005},
        ],
        [
            {"source": "Ancestor", "dest": "Deme2", "rate": 0.01},
            {"source": "Deme2", "dest": "Ancestor", "rate": 0.005},
        ],
        [{"demes": ["Ancestor", "Deme1", "Deme2"], "rate": 1e-5}],
    ],
)
def test_three_deme_sequential_branches(
    second_ancestor, first_start_time, second_start_time, migrations
):
    b = demes.Builder(description="test branch", time_units="generations")
    b.add_deme(name="Ancestor", epochs=[dict(start_size=500, end_time=0)])
    b.add_deme(
        "Deme1",
        epochs=[dict(start_size=100, end_time=0)],
        ancestors=["Ancestor"],
        start_time=first_start_time,
    )
    b.add_deme(
        "Deme2",
        epochs=[dict(start_size=100, end_time=0)],
        ancestors=[second_ancestor],
        start_time=second_start_time,
    )
    for m in migrations:
        b.add_migration(**m)
    resolve_and_run(b)


@run_model_round_trip
class TestMultipleBranches(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(description="test branch",
                               time_units="generations")
        self.b.add_deme(name="Ancestor", epochs=[
                        dict(start_size=1000, end_time=0)])
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
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 1)


@run_model_round_trip
class TestMultipleBranchesWithMigration(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(description="test branch",
                               time_units="generations")
        self.b.add_deme(name="Ancestor", epochs=[
                        dict(start_size=1000, end_time=0)])
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
        self.b.add_migration(source="Deme1", dest="Deme2", rate=1e-6)
        self.b.add_migration(source="Deme2", dest="Deme1", rate=2e-6)
        self.b.add_migration(demes=["Deme1", "Deme3"], rate=0.01)
        self.b.add_migration(source="Deme3", dest="Deme2", rate=1e-6)
        self.b.add_migration(source="Deme2", dest="Deme3", rate=2e-6)
        self.g = self.b.resolve()
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 1)


@run_model_round_trip
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
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 1)


@run_model_round_trip
class TestIslandModel(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(description="island", time_units="generations")
        self.b.add_deme(name="Island1", epochs=[
                        dict(start_size=100, end_time=0)])
        self.b.add_deme(name="Island2", epochs=[
                        dict(start_size=200, end_time=0)])
        self.b.add_migration(source="Island1", dest="Island2", rate=0.01)
        self.b.add_migration(source="Island2", dest="Island1", rate=0.02)
        self.g = self.b.resolve()
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 2)

    def test_demog_attributes(self):
        self.assertTrue(
            self.demog.metadata["burnin_time"]
            == sum(self.demog.initial_sizes_list) * 2
        )


@run_model_round_trip
class TestIslandModelRateChange(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(description="island", time_units="generations")
        self.b.add_deme(name="Island1", epochs=[
                        dict(start_size=100, end_time=0)])
        self.b.add_deme(name="Island2", epochs=[
                        dict(start_size=200, end_time=0)])
        self.b.add_migration(
            source="Island1", dest="Island2", rate=0.01, end_time=500)
        self.b.add_migration(source="Island2", dest="Island1", rate=0.02)
        self.b.add_migration(
            source="Island1", dest="Island2", rate=0.05, start_time=500
        )
        self.g = self.b.resolve()
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 2)

    def test_burnin_time(self):
        self.assertTrue(
            self.demog.metadata["burnin_time"]
            == sum(self.demog.initial_sizes_list) * 2
        )

    def test_total_sim_length(self):
        self.assertEqual(
            self.demog.total_simulation_length,
            self.demog.metadata["burnin_time"] + 500,
        )


@run_model_round_trip
class TestTwoPopMerger(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(
            description="split then merger", time_units="generations"
        )
        self.b.add_deme(name="Ancestral", epochs=[
                        dict(start_size=1000, end_time=1000)])
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
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 1)

    def test_total_sim_length(self):
        self.assertEqual(
            self.demog.total_simulation_length,
            self.demog.metadata["burnin_time"] + 1000,
        )


@run_model_round_trip
class TestFourWayMerger(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(
            description="split then merger", time_units="generations"
        )
        self.b.add_deme(name="Ancestral", epochs=[
                        dict(start_size=1000, end_time=1000)])
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
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 1)

    def test_total_sim_length(self):
        self.assertEqual(
            self.demog.total_simulation_length,
            self.demog.metadata["burnin_time"] + 1000,
        )


@run_model_round_trip
class TestPulseMigration(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.b = demes.Builder(description="pulse", time_units="generations")
        self.b.add_deme(name="deme1", epochs=[
                        dict(start_size=100, end_time=0)])
        self.b.add_deme(name="deme2", epochs=[
                        dict(start_size=100, end_time=0)])
        self.b.add_pulse(sources=["deme1"], dest="deme2",
                         time=100, proportions=[0.2])
        self.g = self.b.resolve()
        self.demog = fwdpy11.discrete_demography.from_demes(self.g, 1)

    def test_total_sim_length(self):
        self.assertEqual(
            self.demog.total_simulation_length,
            self.demog.metadata["burnin_time"] + 100,
        )


def evolve_demes_model(demog) -> fwdpy11.DiploidPopulation:
    pdict = {
        "rates": (0.0, 0.0, 0.0),
        "gvalue": fwdpy11.Multiplicative(2.0),
        "demography": demog,
        "simlen": demog.total_simulation_length,
    }
    initial_sizes = demog.initial_sizes_list
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
    return g


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
    return g


@pytest.mark.parametrize(
    "graph",
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
def test_three_way_continuous_migration_pairwise(graph):
    _ = fwdpy11.discrete_demography.from_demes(graph, 1)


@pytest.mark.parametrize(
    "graph",
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
def test_evolve_three_way_continuous_migration_pairwise(graph):
    demog = fwdpy11.discrete_demography.from_demes(graph, 1)
    pop = evolve_demes_model(demog)
    assert pop.generation == demog.total_simulation_length


def multiple_migrations_delayed():
    b = demes.Builder(description="many migrations", time_units="generations")
    b.add_deme(name="A", epochs=[dict(start_size=100, end_time=0)])
    b.add_deme(name="B", epochs=[dict(start_size=100, end_time=0)])
    b.add_deme(name="C", epochs=[dict(start_size=100, end_time=0)])
    b.add_migration(demes=["A", "B"], rate=0.1)
    b.add_migration(demes=["B", "C"], rate=0.1)
    b.add_migration(demes=["A", "C"], rate=0.1, start_time=100)
    g = b.resolve()
    return g


@pytest.mark.parametrize("graph", [multiple_migrations_delayed()])
def test_multiple_migrations_delayed(graph):
    _ = fwdpy11.discrete_demography.from_demes(graph, 1)


def splits_with_migrations():
    b = demes.Builder(description="splits with migration",
                      time_units="generations")
    b.add_deme(name="A", epochs=[dict(start_size=100, end_time=100)])
    b.add_deme(name="B", ancestors=["A"], epochs=[
               dict(start_size=100, end_time=0)])
    b.add_deme(name="C", ancestors=["A"], epochs=[
               dict(start_size=100, end_time=50)])
    b.add_deme(name="D", ancestors=["C"], epochs=[
               dict(start_size=100, end_time=0)])
    b.add_deme(name="E", ancestors=["C"], epochs=[
               dict(start_size=100, end_time=0)])
    b.add_migration(demes=["B", "C"], rate=0.1)
    b.add_migration(demes=["B", "D"], rate=0.1)
    b.add_migration(demes=["B", "E"], rate=0.1)
    b.add_migration(demes=["D", "E"], rate=0.1)
    g = b.resolve()
    return g


@pytest.mark.parametrize("graph", [splits_with_migrations()])
def test_splits_with_migrations(graph):
    _ = fwdpy11.discrete_demography.from_demes(graph, 1)


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
    _ = fwdpy11.discrete_demography.from_demes(g, 1)


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
                    self.sizes[key] = [DemeSizeAtTime(
                        when=pop.generation, size=value)]
                else:
                    self.sizes[key].append(
                        DemeSizeAtTime(when=pop.generation, size=value)
                    )

    pdict = {
        "gvalue": fwdpy11.Multiplicative(2.0),
        "rates": (0, 0, 0),
        "demography": model,
        "simlen": model.total_simulation_length,
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
            i for i in range(111, model.total_simulation_length + 1)
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

    pdict["simlen"] = model.total_simulation_length - when
    params = fwdpy11.ModelParams(**pdict)

    fwdpy11.evolvets(rng, pop, params, 100)

    counts = np.unique(np.array(pop.diploid_metadata)
                       ["deme"], return_counts=True)
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
    pdict2["simlen"] = model.total_simulation_length - when
    params2 = fwdpy11.ModelParams(**pdict2)

    fwdpy11.evolvets(rng, pop, params2, 100)

    counts = np.unique(np.array(pop.diploid_metadata)
                       ["deme"], return_counts=True)
    assert counts[1][0] == 500, f"{counts}"
    assert counts[1][1] == 200, f"{counts}"


@pytest.mark.parametrize("burnin", [1, 0])
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
        "simlen": model.total_simulation_length,
    }
    params = fwdpy11.ModelParams(**pdict)
    initial_sizes = model.initial_sizes_list
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
                    self.sizes[key] = [DemeSizeAtTime(
                        when=pop.generation, size=value)]
                else:
                    self.sizes[key].append(
                        DemeSizeAtTime(when=pop.generation, size=value)
                    )

    pdict = {
        "gvalue": fwdpy11.Multiplicative(2.0),
        "rates": (0, 0, 0),
        "demography": model,
        "simlen": model.total_simulation_length,
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
                    self.sizes[key] = [DemeSizeAtTime(
                        when=pop.generation, size=value)]
                else:
                    self.sizes[key].append(
                        DemeSizeAtTime(when=pop.generation, size=value)
                    )

    pdict = {
        "gvalue": fwdpy11.Multiplicative(2.0),
        "rates": (0, 0, 0),
        "demography": model,
        "simlen": model.total_simulation_length,
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
                TrackedInfo(pop.generation, deme,
                            selfed_by_deme[deme], size, pop.N)
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
        - {end_time: 1, selfing_rate: 0.0, start_size: 100}
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
        "simlen": demog.total_simulation_length,
    }
    params = fwdpy11.ModelParams(**pdict)

    initial_sizes = demog.initial_sizes_list
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

    initial_sizes = demog.initial_sizes_list
    pop = fwdpy11.DiploidPopulation(initial_sizes, 1.0)
    recorder(pop, None)
    rng = fwdpy11.GSLrng(90210)
    fwdpy11.evolvets(rng, pop, params, 100, recorder=recorder)

    pdict["simlen"] = demog.total_simulation_length - 1
    fwdpy11.evolvets(
        rng, pop, params, 100, recorder=recorder,
    )

    validate(recorder)


@pytest.fixture
def multiple_pulse_source_setup():
    b = demes.Builder(description="pulse", time_units="generations")
    b.add_deme(name="deme0", epochs=[dict(start_size=100, end_time=0)])
    b.add_deme(name="deme1", epochs=[dict(start_size=100, end_time=0)])
    b.add_deme(name="deme2", epochs=[dict(start_size=100, end_time=0)])
    b.add_pulse(
        sources=["deme0", "deme1"], dest="deme2", time=100, proportions=[0.2, 0.25]
    )
    g = b.resolve()
    demog = fwdpy11.discrete_demography.from_demes(g, 10)

    return demog


def test_multiple_pulse_source(multiple_pulse_source_setup):
    demog = multiple_pulse_source_setup
    assert (
        demog.total_simulation_length == demog.metadata["burnin_time"] + 100
    )


@pytest.fixture
def two_demes_migration_rate_changes_setup():
    b = demes.Builder(time_units="generations")
    b.add_deme(name="A", epochs=[dict(start_size=100)])
    b.add_deme(name="B", epochs=[dict(start_size=100)])
    b.add_migration(demes=["A", "B"], rate=0.0, end_time=30)
    b.add_migration(demes=["A", "B"], rate=0.5, start_time=30, end_time=20)
    b.add_migration(demes=["A", "B"], rate=0.0, start_time=20, end_time=10)
    b.add_migration(demes=["A", "B"], rate=0.5, start_time=10)
    g = b.resolve()
    return g


def test_two_demes_migration_rate_changes(two_demes_migration_rate_changes_setup):
    # test to check that the last 10 generation have the high migration rates
    # not 9 generations, and not 11 generations

    @dataclass
    class ParentDemesAtTime:
        when: int
        parents: list

    class ParentDemes(object):
        def __init__(self):
            self.parents = dict()

        def __call__(self, pop, _):
            deme_indexes, deme_sizes = pop.deme_sizes()
            for deme_idx in deme_indexes:
                if deme_idx not in self.parents:
                    self.parents[deme_idx] = [
                        ParentDemesAtTime(when=pop.generation, parents=[0, 0])
                    ]
                else:
                    self.parents[deme_idx].append(
                        ParentDemesAtTime(when=pop.generation, parents=[0, 0])
                    )
            for metadata in pop.diploid_metadata:
                parents = metadata.parents
                for parent in parents:
                    self.parents[metadata.deme][-1].parents[
                        pop.diploid_metadata[parent].deme
                    ] += 1

    demog = fwdpy11.discrete_demography.from_demes(
        two_demes_migration_rate_changes_setup, 1
    )

    pdict = {
        "gvalue": fwdpy11.Multiplicative(2.0),
        "rates": (0, 0, 0),
        "demography": demog,
        "simlen": demog.total_simulation_length,
    }
    params = fwdpy11.ModelParams(**pdict)
    initial_sizes = demog.initial_sizes_list
    pop = fwdpy11.DiploidPopulation([initial_sizes[0], initial_sizes[1]], 1.0)
    rng = fwdpy11.GSLrng(90210)
    recorder = ParentDemes()
    fwdpy11.evolvets(rng, pop, params, 100, recorder=recorder)
    for ii in range(1, 11):
        assert (
            recorder.parents[0][-ii].parents[0] > 0
            and recorder.parents[0][-ii].parents[1] > 0
        )
        assert (
            recorder.parents[0][-ii].parents[0] > 0
            and recorder.parents[0][-ii].parents[1] > 0
        )
        assert (
            recorder.parents[0][-ii - 20].parents[0] > 0
            and recorder.parents[0][-ii - 20].parents[1] > 0
        )
        assert (
            recorder.parents[0][-ii - 20].parents[0] > 0
            and recorder.parents[0][-ii - 20].parents[1] > 0
        )
        assert recorder.parents[0][-ii - 10].parents[0] == 200
        assert recorder.parents[1][-ii - 10].parents[1] == 200
        assert recorder.parents[0][-ii - 30].parents[0] == 200
        assert recorder.parents[1][-ii - 30].parents[1] == 200


def test_split_with_existence_n_way_migration():
    b = demes.Builder(time_units="generations")
    b.add_deme(name="anc", epochs=[dict(start_size=100, end_time=200)])
    b.add_deme(name="A", ancestors=["anc"], epochs=[dict(start_size=100)])
    b.add_deme(name="B", ancestors=["anc"], epochs=[
               dict(start_size=100, end_time=100)])
    b.add_deme(name="C", ancestors=["B"], epochs=[dict(start_size=100)])
    b.add_deme(name="D", ancestors=["B"], epochs=[dict(start_size=100)])
    b.add_deme(name="E", ancestors=["B"], epochs=[dict(start_size=100)])
    b.add_migration(demes=["A", "C", "D", "E"], rate=0.1)
    resolve_and_run(b, burnin=1)


def test_non_integer_size_change():
    b = demes.Builder()
    b.add_deme("A", epochs=[{"start_size": 100}])
    b.add_deme("B", ancestors=["A"], start_time=10,
               epochs=[{"start_size": 100.5}])
    g = b.resolve()

    _ = fwdpy11.discrete_demography.from_demes(g, round_non_integer_sizes=True)


def test_non_integer_initial_epoch_size():
    b = demes.Builder()
    b.add_deme("A", epochs=[{"start_size": 100.5}])
    g = b.resolve()

    demog = fwdpy11.discrete_demography.from_demes(
        g, round_non_integer_sizes=True)

    # initialize the ancestral population
    initial_sizes = demog.initial_sizes_list
    _ = fwdpy11.DiploidPopulation(
        initial_sizes, 1)


def test_very_short_epoch():
    b = demes.Builder()
    b.add_deme("A", epochs=[{"start_size": 100}])
    b.add_deme(
        "B",
        ancestors=["A"],
        start_time=10,
        epochs=[
            {"start_size": 100, "end_time": 1e-3},
            {"start_size": 200, "end_time": 0.0},
        ],
    )
    g = b.resolve()

    with pytest.raises(ValueError):
        _ = fwdpy11.discrete_demography.from_demes(g)


def test_epoch_rounding_01():
    yaml = """
time_units: generations
demes:
- name: bad
  epochs:
  - {end_time: 1.75, start_size: 1}
  - {end_time: 1.25, start_size: 2}
  - {end_time: 0.5, start_size: 3}
  - {end_time: 0, start_size: 4}
"""
    graph = demes.loads(yaml)
    with pytest.raises(ValueError):
        _ = fwdpy11.discrete_demography.from_demes(graph, burnin=0)


def test_epoch_rounding_02():
    yaml = """
time_units: generations
demes:
- name: bad
  epochs:
  - {end_time: 1.5, start_size: 1}
  - {end_time: 0.4, start_size: 2}
  - {end_time: 0, start_size: 3}
"""
    graph = demes.loads(yaml)
    with pytest.raises(ValueError):
        _ = fwdpy11.discrete_demography.from_demes(graph, burnin=0)


def test_ambiguous_pulses():
    yaml = """
time_units: generations
demes:
 - name: A
   epochs:
    - start_size: 50
 - name: B
   epochs:
    - start_size: 50
pulses:
 - sources: [A]
   dest: B
   proportions: [0.9]
   time: 25
 - sources: [A]
   dest: B
   proportions: [0.9]
   time: 25
"""
    with pytest.warns(UserWarning):
        graph = demes.loads(yaml)
    with pytest.warns(UserWarning):
        _ = fwdpy11.discrete_demography.from_demes(graph, burnin=0)


def test_deme_sort_order():
    yaml = """
description: single deme model
time_units: generations
demes:
 - name: B
   epochs:
    - start_size: 100
 - name: A
   epochs:
    - start_size: 100
"""
    graph = demes.loads(yaml)
    model = fwdpy11.discrete_demography.from_demes(graph, burnin=0)
    assert model.deme_labels[0] == "B"
    assert model.deme_labels[1] == "A"


def test_residual_selfing_exception():
    yaml = """
description: single deme model
time_units: generations
demes:
 - name: B
   epochs:
    - start_size: 1
"""
    graph = demes.loads(yaml)
    model = fwdpy11.discrete_demography.from_demes(graph, burnin=10)
    pop = fwdpy11.DiploidPopulation([1], 1.0)
    pdict = {"recregions": [],
             "nregions": [],
             "rates": (0, 0, 0),
             "gvalue": fwdpy11.Multiplicative(2.),
             "simlen": 10,
             "demography": model,
             "allow_residual_selfing": False
             }
    params = fwdpy11.ModelParams(**pdict)
    rng = fwdpy11.GSLrng(42)
    with pytest.raises(fwdpy11.DemographyError):
        with pytest.warns(UserWarning):
            fwdpy11.evolvets(rng, pop, params, 5)
