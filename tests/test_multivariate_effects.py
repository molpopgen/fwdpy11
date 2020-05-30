import unittest

import numpy as np

import fwdpy11
import poptools
import testMultivariateGSSmo


def set_up_quant_trait_model():
    N = 1000
    rho = 1000.0
    r = rho / (4 * N)

    timepoints = np.array([0], dtype=np.uint32)
    ntraits = 4
    optima = np.array(np.zeros(len(timepoints) * ntraits))
    optima = optima.reshape((len(timepoints), ntraits))
    PO = fwdpy11.PleiotropicOptima
    po = []
    for i, j in enumerate(timepoints):
        po.append(PO(when=int(j), optima=optima[i, :], VS=1.0))
    GSSmo = fwdpy11.MultivariateGSSmo(po)
    cmat = np.identity(ntraits)
    np.fill_diagonal(cmat, 0.1)
    a = fwdpy11.StrictAdditiveMultivariateEffects(ntraits, 0, GSSmo)
    p = {
        "nregions": [],
        "sregions": [fwdpy11.MultivariateGaussianEffects(0, 1, 1, cmat)],
        "recregions": [fwdpy11.Region(0, 1, 1)],
        "rates": (0.0, 0.001, r),
        "gvalue": a,
        "prune_selected": False,
        "demography": fwdpy11.DiscreteDemography(),
        "simlen": 10 * N,
    }
    params = fwdpy11.ModelParams(**p)
    rng = fwdpy11.GSLrng(101 * 45 * 110 * 210)
    pop = fwdpy11.DiploidPopulation(N, 1.0)
    return params, rng, pop, ntraits


class TestMultivariateGSSmo(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.params, self.rng, self.pop, self.ntraits = set_up_quant_trait_model()

    def test_github_issue_310(self):
        fwdpy11.evolvets(self.rng, self.pop, self.params, 100)
        nonmutants = 0
        md = np.array(self.pop.diploid_metadata, copy=False)
        for i, j in enumerate(self.pop.diploids):
            n0 = len(self.pop.haploid_genomes[j.first].smutations)
            n1 = len(self.pop.haploid_genomes[j.second].smutations)
            if md["w"][i] != 1.0:
                self.assertTrue(n0 + n1 > 0)
            elif n1 + n0 == 0:
                nonmutants += 1

            if nonmutants == 0:
                self.fail("Test invalid: zero individuals were mutation-free")


class TestMultivariateGSSmoCPP(unittest.TestCase):
    """
    Tests of C++ API
    """

    def setUp(self):
        self.pop = fwdpy11.DiploidPopulation(100, 1.0)
        self.ntraits = 4

    def test_multiple_optima(self):
        timepoints = np.array([0, 5 * self.pop.N])
        optima = np.array(np.zeros(len(timepoints) * self.ntraits))
        optima = optima.reshape((len(timepoints), self.ntraits))
        optima[1, 0] = 1
        PO = fwdpy11.PleiotropicOptima
        po = []
        for i, t in enumerate(timepoints):
            po.append(PO(when=int(t), optima=optima[i, :], VS=1.0))
        GSSmo = fwdpy11.MultivariateGSSmo(po)

        co = testMultivariateGSSmo.get_optima(GSSmo)
        self.assertTrue(all([i == 0 for i in co]))
        poptools.change_generation(self.pop, timepoints[1])
        testMultivariateGSSmo.update(self.pop, GSSmo)
        co = testMultivariateGSSmo.get_optima(GSSmo)
        self.assertEqual(co[0], 1)
        self.assertTrue(all([i == 0 for i in co[1:]]))

        poptools.change_generation(self.pop, timepoints[1] + 10)
        testMultivariateGSSmo.update(self.pop, GSSmo)
        co2 = testMultivariateGSSmo.get_optima(GSSmo)

        self.assertTrue(co == co2)


if __name__ == "__main__":
    unittest.main()
