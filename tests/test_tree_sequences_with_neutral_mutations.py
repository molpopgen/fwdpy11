import unittest
import fwdpy11
import numpy as np
import copy
from test_tree_sequences import set_up_quant_trait_model
from test_tree_sequences import set_up_standard_pop_gen_model


def _count_mutations_from_diploids(pop):
    mc = np.zeros(len(pop.mutations), dtype=np.uint32)
    for d in pop.diploids:
        for k in pop.haploid_genomes[d.first].mutations:
            mc[k] += 1
        for k in pop.haploid_genomes[d.first].smutations:
            mc[k] += 1
        for k in pop.haploid_genomes[d.second].mutations:
            mc[k] += 1
        for k in pop.haploid_genomes[d.second].smutations:
            mc[k] += 1
    return mc


class TestKeepFixations(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.params, self.rng, self.pop = set_up_quant_trait_model()
        self.params.rates = (1e-3, self.params.rates[1], self.params.rates[2])
        self.params.nregions = [fwdpy11.Region(0, 1, 1)]

    def test_mutation_counts_with_indexing(self):
        """
        Tests atypical use case where neutral mutations are placed into genomes
        """
        fwdpy11.evolvets(self.rng, self.pop, self.params, 100,
                         put_neutral_variants_in_genomes=True)
        mc = _count_mutations_from_diploids(self.pop)
        self.assertTrue(np.array_equal(
            mc, np.array(self.pop.mcounts, dtype=np.uint32)))

    def test_mutation_counts_with_indexing_suppressed(self):
        """
        Tests atypical use case where neutral mutations are placed into genomes
        """
        fwdpy11.evolvets(self.rng, self.pop, self.params, 100,
                         suppress_table_indexing=True,
                         put_neutral_variants_in_genomes=True)
        mc = _count_mutations_from_diploids(self.pop)
        self.assertTrue(np.array_equal(
            mc, np.array(self.pop.mcounts, dtype=np.uint32)))

    def test_mutation_counts_with_indexing_suppressed_no_neutral_muts_in_genomes(self):
        """
        A sim w/ and w/o putting neutral variants in genomes
        should give the same mutation counts.
        """
        pop2 = copy.deepcopy(self.pop)
        rng = fwdpy11.GSLrng(101*45*110*210)  # Use same seed!!!
        self.params.prune_selected = False
        params = copy.deepcopy(self.params)
        fwdpy11.evolvets(rng, pop2, params, 100,
                         put_neutral_variants_in_genomes=True,
                         suppress_table_indexing=True)
        fwdpy11.evolvets(self.rng, self.pop, self.params, 100,
                         put_neutral_variants_in_genomes=False,
                         suppress_table_indexing=True)
        ti = fwdpy11.TreeIterator(
            self.pop.tables, [i for i in range(2*self.pop.N)])
        mc = _count_mutations_from_diploids(self.pop)
        for t in ti:
            for m in t.mutations():
                # Have to skip neutral mutations b/c they won't
                # end up in mc b/c it is obtained from genomes
                if pop2.mutations[m.key].neutral is False:
                    self.assertEqual(mc[m.key], self.pop.mcounts[m.key])
                self.assertEqual(t.leaf_counts(m.node), pop2.mcounts[m.key])


class TestPruneFixations(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.params, self.rng, self.pop = set_up_standard_pop_gen_model()
        self.params.rates = (1e-3, self.params.rates[1], self.params.rates[2])
        self.params.nregions = [fwdpy11.Region(0, 1, 1)]

    def test_mutation_counts_with_indexing(self):
        """
        Tests atypical use case where neutral mutations are placed into genomes
        """
        fwdpy11.evolvets(self.rng, self.pop, self.params, 100,
                         put_neutral_variants_in_genomes=True)
        mc = _count_mutations_from_diploids(self.pop)
        self.assertTrue(np.array_equal(
            mc, np.array(self.pop.mcounts, dtype=np.uint32)))

    def test_mutation_counts_with_indexing_suppressed(self):
        """
        Tests atypical use case where neutral mutations are placed into genomes
        """
        fwdpy11.evolvets(self.rng, self.pop, self.params, 100,
                         suppress_table_indexing=True,
                         put_neutral_variants_in_genomes=True)
        mc = _count_mutations_from_diploids(self.pop)
        self.assertTrue(np.array_equal(
            mc, np.array(self.pop.mcounts, dtype=np.uint32)))

    def test_mutation_counts_with_indexing_suppressed_no_neutral_muts_in_genomes(self):
        """
        A sim w/ and w/o putting neutral variants in genomes
        should give the same mutation counts.
        """
        pop2 = copy.deepcopy(self.pop)
        rng = fwdpy11.GSLrng(666**2)  # Use same seed!!!
        params = copy.deepcopy(self.params)
        fwdpy11.evolvets(rng, pop2, params, 100,
                         put_neutral_variants_in_genomes=True,
                         suppress_table_indexing=True)
        fwdpy11.evolvets(self.rng, self.pop, self.params, 100,
                         put_neutral_variants_in_genomes=False,
                         suppress_table_indexing=True)
        ti = fwdpy11.TreeIterator(
            self.pop.tables, [i for i in range(2*self.pop.N)])
        mc = _count_mutations_from_diploids(self.pop)
        for t in ti:
            for m in t.mutations():
                # Have to skip neutral mutations b/c they won't
                # end up in mc b/c it is obtained from genomes
                if pop2.mutations[m.key].neutral is False:
                    self.assertEqual(mc[m.key], self.pop.mcounts[m.key])
                self.assertEqual(t.leaf_counts(m.node), pop2.mcounts[m.key])


if __name__ == "__main__":
    unittest.main()
