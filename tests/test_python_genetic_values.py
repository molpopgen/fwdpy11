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
import unittest
import fwdpy11 as fp11
import fwdpy11.fitness as fp11w
from quick_pops import quick_nonneutral_slocus
import fwdpy11.python_genetic_values as pgv


def additive(dip, gametes, mutations):
    """
    Pure Python implementation of additive fitness.
    """
    rv = 0.0
    for i in gametes[dip.first].smutations:
        rv += mutations[i].s
    for i in gametes[dip.second].smutations:
        rv += mutations[i].s
    return max(0.0, 1 + rv)


class testCustomAdditive(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.pop = quick_nonneutral_slocus()
        self.acpp = fp11w.SlocusAdditive(2.0)
        self.a = pgv.GeneticValue(additive)

    def testDiploidFitnesses(self):
        """
        Compare our Python implementation to
        our C++ implementation.
        """
        for i in self.pop.diploids:
            self.assertEqual(self.acpp(i, self.pop), self.a(i, self.pop))


class testAdditiveEvolution(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        from fwdpy11.ezparams import mslike
        from fwdpy11.regions import ExpS
        from fwdpy11.model_params import SlocusParams
        self.pop = fp11.SlocusPop(1000)
        self.pdict = mslike(
            self.pop, simlen=100, dfe=ExpS(0, 1, 1, -0.1),
            pneutral=0.95)
        self.pdict['gvalue'] = pgv.GeneticValue(additive)
        self.params = SlocusParams(**self.pdict)
        self.rng = fp11.GSLrng(42)

    def testValidate(self):
        try:
            self.params.validate()
        except:
            self.fail("Parameters failed to validate")

    def testEvolve(self):
        from fwdpy11.wright_fisher import evolve
        try:
            evolve(self.rng, self.pop, self.params)
        except:
            self.fail("unexpected exception")


if __name__ == "__main__":
    unittest.main()
