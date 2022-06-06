import fwdpy11


def test_basic_semantics():
    g = fwdpy11._fwdpy11._ll_ForwardDemesGraph()
    _ = g.add_deme("bananas", 0, 0)
