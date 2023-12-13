import demes

import fwdpy11


def into_graph(yaml, burnin=100, burnin_is_exact=True) -> fwdpy11.ForwardDemesGraph:
    graph = demes.load(yaml)
    return fwdpy11.ForwardDemesGraph.from_demes(
        graph, burnin=burnin, burnin_is_exact=burnin_is_exact
    )


# This function is circular logic as we are repeating
# the current implementation of to_forwards_time
def validate_deme_times(demog):
    for deme in demog.graph.demes:
        if deme.start_time != float("inf"):
            assert (
                demog.to_forwards_time(deme.start_time)
                == demog.final_generation - deme.start_time
            )
        else:
            assert demog.to_forwards_time(deme.start_time) is None
        assert (
            demog.to_forwards_time(deme.end_time)
            == demog.final_generation - deme.end_time
        )
    deme_times = demog.deme_time_intervals()
    assert max([i.end_time for i in deme_times]) == demog.final_generation + 1


def test_single_pulse():
    yaml = "tests/demes_event_examples/single_pulse.yaml"
    demog = into_graph(yaml)
    assert demog.to_forwards_time(demog.final_generation) == 0
    assert demog.to_forwards_time(-1) is None
    assert demog.to_forwards_time(float("inf")) is None
    assert demog.to_forwards_time(10) == 100
    assert demog.to_backwards_time(10) == 100
    assert demog.to_backwards_time(0) == demog.final_generation
    validate_deme_times(demog)


def test_burst_of_migration():
    yaml = "tests/demes_event_examples/burst_of_migration.yaml"
    demog = into_graph(yaml)
    assert demog.to_forwards_time(demog.final_generation) == 0
    assert demog.to_forwards_time(-1) is None
    assert demog.to_forwards_time(10) == 100
    assert demog.to_backwards_time(0) == demog.final_generation
    validate_deme_times(demog)


def test_deme_existence():
    yaml = "tests/demes_event_examples/deme_existence.yaml"
    demog = into_graph(yaml)
    assert demog.to_forwards_time(demog.final_generation) == 0
    assert demog.to_forwards_time(-1) is None
    assert demog.to_forwards_time(10) == 140
    assert demog.to_backwards_time(0) == demog.final_generation
    validate_deme_times(demog)
    assert demog.to_forwards_time(demog.graph.demes[0].start_time) is None
    assert demog.to_forwards_time(demog.graph.demes[1].start_time) == 100
    assert demog.to_forwards_time(demog.graph.demes[2].start_time) == 125
    assert demog.to_forwards_time(demog.graph.demes[0].end_time) == 100
    assert demog.to_forwards_time(demog.graph.demes[1].end_time) == 125
    assert demog.to_forwards_time(demog.graph.demes[2].end_time) == 150
    assert (
        demog.to_forwards_time(demog.graph.demes[2].end_time) == demog.final_generation
    )


def test_model_that_ends_before_zero():
    yaml = """
    time_units: generations
    demes:
     - name: deme0
       epochs:
        - start_size: 100
          end_time: 50
        - start_size: 100
          end_time: 3
    """
    graph = demes.loads(yaml)
    demog = fwdpy11.ForwardDemesGraph.from_demes(
        graph, burnin=100, burnin_is_exact=True
    )
    assert demog.final_generation == 147
    assert demog.to_forwards_time(3) == demog.final_generation
    assert demog.to_backwards_time(demog.final_generation) == 3
    assert demog.to_backwards_time(0) == 150
    assert demog.to_backwards_time(3) == 147
