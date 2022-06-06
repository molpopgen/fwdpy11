import fwdpy11


# FIXME: remove this test
# when we have a working public API
def test_basic_semantics():
    g = fwdpy11._fwdpy11._ll_ForwardDemesGraph()
    d = g._add_deme("bananas", 0, 0)
    d._add_epoch(
        end_time=0,
        start_size=100,
        end_size=100,
        cloning_rate=0.0,
        selfing=fwdpy11._fwdpy11._wright_fisher_selfing(),
        size_function=fwdpy11._fwdpy11._constant_size_function(),
    )
