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

# TYPETEST

import fwdpy11
import unittest


class testSlocusPopCreate(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.mutations = fwdpy11.VecMutation()
        self.gametes = fwdpy11.VecGamete()
        self.diploids = fwdpy11.VecDiploid()
        self.mutations.append(fwdpy11.Mutation(0.1, -0.01, 1.0, 0, 0))
        self.gametes.append(fwdpy11.Gamete(
            (2, fwdpy11.VecUint32([]), fwdpy11.VecUint32([0]))))
        self.diploids.append(fwdpy11.SingleLocusDiploid(0, 0))

    def testConstruction(self):
        pop = fwdpy11.SlocusPop(self.diploids, self.gametes, self.mutations)
        self.assertTrue(pop.N, 1)

    def testStaticMethod(self):
        pop = fwdpy11.SlocusPop.create(
            self.diploids, self.gametes, self.mutations)
        self.assertTrue(type(pop) is fwdpy11.SlocusPop)

    def testStaticMethodWithFixations(self):
        ftimes = fwdpy11.VecUint32([1])
        pop = fwdpy11.SlocusPop.create(self.diploids,
                self.gametes,self.mutations,self.mutations,ftimes,2)


if __name__ == "__main__":
    unittest.main()
