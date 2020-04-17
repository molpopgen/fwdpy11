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

import fwdpy11


class testDiploidPopulationAddMutations(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.pop = fwdpy11.DiploidPopulation(1000)
        self.mvec = fwdpy11.MutationVector()
        self.nmut = fwdpy11.Mutation(0.1, 0.0, 0.0, 0, 0)
        self.smut = fwdpy11.Mutation(1.2, -0.01, 0.0, 0, 0)
        self.smut_vec = fwdpy11.Mutation(1.3, 0.0, 0.0, 0, [-1.0, 0], [1, 1], 0)

    def testAddOneNeutralMutation(self):
        self.mvec.append(self.nmut)
        m = self.pop.add_mutations(self.mvec, [0], [0])
        self.assertEqual(len(m), 1)
        for i in m:
            self.assertEqual(self.pop.mutations[i].neutral, True)
            self.assertEqual(self.pop.mcounts[i], 1)
        self.assertEqual(
            len(self.pop.haploid_genomes[self.pop.diploids[0].first].mutations), 1
        )
        self.assertEqual(
            len(self.pop.haploid_genomes[self.pop.diploids[0].first].smutations), 0
        )
        self.assertEqual(
            len(self.pop.haploid_genomes[self.pop.diploids[0].second].mutations), 0
        )
        self.assertEqual(
            len(self.pop.haploid_genomes[self.pop.diploids[0].second].smutations), 0
        )

    def testAddOneSelectedMutation(self):
        self.mvec.append(self.smut)
        m = self.pop.add_mutations(self.mvec, [0], [2])
        self.assertEqual(len(m), 1)
        for i in m:
            self.assertEqual(self.pop.mutations[i].neutral, False)
            self.assertEqual(self.pop.mcounts[i], 2)
        self.assertEqual(
            len(self.pop.haploid_genomes[self.pop.diploids[0].first].mutations), 0
        )
        self.assertEqual(
            len(self.pop.haploid_genomes[self.pop.diploids[0].first].smutations), 1
        )
        self.assertEqual(
            len(self.pop.haploid_genomes[self.pop.diploids[0].second].mutations), 0
        )
        self.assertEqual(
            len(self.pop.haploid_genomes[self.pop.diploids[0].second].smutations), 1
        )

    def testAddOneSelectedMutationWithVecEffects(self):
        self.mvec.append(self.smut_vec)
        m = self.pop.add_mutations(self.mvec, [0], [2])

        # The addition of mutations to populations works
        # via move construction, meaning the vector contents
        # of the input must be empty:
        self.assertEqual(len(self.mvec[0].esizes), 0)
        self.assertEqual(len(self.mvec[0].heffects), 0)

        self.assertEqual(len(m), 1)
        for i in m:
            self.assertEqual(self.pop.mutations[i].neutral, False)
            self.assertEqual(self.pop.mcounts[i], 2)
        self.assertEqual(
            len(self.pop.haploid_genomes[self.pop.diploids[0].first].mutations), 0
        )
        self.assertEqual(
            len(self.pop.haploid_genomes[self.pop.diploids[0].first].smutations), 1
        )
        self.assertEqual(
            len(self.pop.haploid_genomes[self.pop.diploids[0].second].mutations), 0
        )
        self.assertEqual(
            len(self.pop.haploid_genomes[self.pop.diploids[0].second].smutations), 1
        )

    def testAddMultipleMutationsSimple(self):
        self.mvec.append(self.smut)
        self.mvec.append(self.nmut)
        # Add mut to every other individual,
        # and they are all homozygous
        ind = [i for i in range(0, self.pop.N, 2)]
        m = self.pop.add_mutations(self.mvec, ind, [2] * len(ind))
        for i in m:
            self.assertEqual(self.pop.mcounts[i], 2 * len(ind))

    def testAddMultipleMutationsComplex(self):
        self.mvec.append(self.smut)
        self.mvec.append(self.nmut)
        ind1 = [i for i in range(0, int(self.pop.N / 2), 2)]
        ind2 = [i for i in range(int(self.pop.N / 2), self.pop.N, 3)]
        c1 = [0] * len(ind1)
        c2 = [1] * len(ind2)

        x = self.pop.add_mutations(self.mvec, ind1 + ind2, c1 + c2)
        self.assertEqual(len(x), 2)
        self.assertEqual(len(self.pop.mutations), 2)
        for i in x:
            self.assertEqual(self.pop.mcounts[i], len(ind1) + len(ind2))
        for i in ind1:
            self.assertEqual(
                len(self.pop.haploid_genomes[self.pop.diploids[i].first].mutations), 1
            )
            self.assertEqual(
                len(self.pop.haploid_genomes[self.pop.diploids[i].first].smutations), 1
            )
            self.assertEqual(
                len(self.pop.haploid_genomes[self.pop.diploids[i].second].mutations), 0
            )
            self.assertEqual(
                len(self.pop.haploid_genomes[self.pop.diploids[i].second].smutations), 0
            )
        for i in ind2:
            self.assertEqual(
                len(self.pop.haploid_genomes[self.pop.diploids[i].first].mutations), 0
            )
            self.assertEqual(
                len(self.pop.haploid_genomes[self.pop.diploids[i].first].smutations), 0
            )
            self.assertEqual(
                len(self.pop.haploid_genomes[self.pop.diploids[i].second].mutations), 1
            )
            self.assertEqual(
                len(self.pop.haploid_genomes[self.pop.diploids[i].second].smutations), 1
            )

    def testAddSameMutationPosition(self):
        self.mvec.append(self.nmut)
        m = self.pop.add_mutations(self.mvec, [0], [0])
        self.mvec = fwdpy11.MutationVector()
        self.mvec.append(self.nmut)
        with self.assertRaises(ValueError):
            m = self.pop.add_mutations(self.mvec, [0], [0])

    # Test errors that we expect from fwdpp.
    # Note: these should all be ValueError,
    # as they reflect bad input.  We should fix
    # that upstream.

    def testUnEqualListLengths(self):
        self.mvec.append(self.nmut)
        with self.assertRaises(ValueError):
            m = self.pop.add_mutations(self.mvec, [0], [0, 1])

    def testInvalidIndividual(self):
        self.mvec.append(self.nmut)
        with self.assertRaises(IndexError):
            m = self.pop.add_mutations(self.mvec, [self.pop.N], [0])

    def testInvalidHaploidGenome(self):
        self.mvec.append(self.nmut)
        with self.assertRaises(IndexError):
            m = self.pop.add_mutations(self.mvec, [53], [3])


if __name__ == "__main__":
    unittest.main()
