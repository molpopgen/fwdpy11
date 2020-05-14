import copy
import unittest

import numpy as np

import fwdpy11
from test_tree_sequences import set_up_quant_trait_model, set_up_standard_pop_gen_model


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


def _compare_counts_for_nonneutral_variants(pop, mc):
    w = np.array(
        [i for i, j in enumerate(pop.mutations) if j.neutral is False], dtype=np.uint32
    )
    return all(pop.mcounts[w] == mc[w]) is True


class TestKeepFixations(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.params, self.rng, self.pop = set_up_quant_trait_model()
        pdict = self.params.asdict()
        pdict["rates"] = (
            1e-3,
            self.params.rates.selected_mutation_rate,
            self.params.rates.recombination_rate,
        )
        pdict["nregions"] = [fwdpy11.Region(0, 1, 1)]
        self.params = fwdpy11.ModelParams(**pdict)

    def test_mutation_counts_with_indexing(self):
        """
        Tests atypical use case where neutral mutations are placed into genomes
        """
        fwdpy11.evolvets(self.rng, self.pop, self.params, 100)
        mc = _count_mutations_from_diploids(self.pop)
        self.assertTrue(_compare_counts_for_nonneutral_variants(self.pop, mc))

    def test_mutation_counts_with_indexing_suppressed(self):
        """
        Tests atypical use case where neutral mutations are placed into genomes
        """
        fwdpy11.evolvets(
            self.rng, self.pop, self.params, 100, suppress_table_indexing=True
        )
        mc = _count_mutations_from_diploids(self.pop)
        self.assertTrue(_compare_counts_for_nonneutral_variants(self.pop, mc))

    def test_mutation_counts_with_indexing_suppressed_no_neutral_muts_in_genomes(self):
        """
        A sim w/ and w/o putting neutral variants in genomes
        should give the same mutation counts.
        """
        pop2 = copy.deepcopy(self.pop)
        rng = fwdpy11.GSLrng(101 * 45 * 110 * 210)  # Use same seed!!!
        pdict = self.params.asdict()
        pdict["prune_selected"] = False
        self.params = fwdpy11.ModelParams(**pdict)
        params = copy.deepcopy(self.params)
        fwdpy11.evolvets(rng, pop2, params, 100, suppress_table_indexing=True)
        fwdpy11.evolvets(
            self.rng, self.pop, self.params, 100, suppress_table_indexing=True
        )
        ti = fwdpy11.TreeIterator(self.pop.tables, [i for i in range(2 * self.pop.N)])
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
        pdict = self.params.asdict()
        pdict["rates"] = (
            1e-3,
            self.params.rates.selected_mutation_rate,
            self.params.rates.recombination_rate,
        )
        pdict["nregions"] = [fwdpy11.Region(0, 1, 1)]
        self.params = fwdpy11.ModelParams(**pdict)

    @unittest.skip("doesn't test anything...")
    def test_no_mutations(self):
        pdict = self.params.asdict()
        pdict["rates"] = (0, 0, 0)
        self.params = fwdpy11.ModelParams(**pdict)
        fwdpy11.evolvets(self.rng, self.pop, self.params, 100)

    def test_mutation_counts_with_indexing(self):
        fwdpy11.evolvets(self.rng, self.pop, self.params, 100)
        mc = _count_mutations_from_diploids(self.pop)
        self.assertTrue(_compare_counts_for_nonneutral_variants(self.pop, mc))
        self.assertEqual(len(np.where(mc == 2 * self.pop.N)[0]), 0)

    def test_mutation_counts_with_indexing_suppressed(self):
        fwdpy11.evolvets(
            self.rng, self.pop, self.params, 100, suppress_table_indexing=True
        )
        mc = _count_mutations_from_diploids(self.pop)
        self.assertTrue(_compare_counts_for_nonneutral_variants(self.pop, mc))
        self.assertEqual(len(np.where(mc == 2 * self.pop.N)[0]), 0)


class TestNeutralMutRegions(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.params, self.rng, self.pop = set_up_standard_pop_gen_model()
        pdict = self.params.asdict()
        pdict["rates"] = (
            1e-3,
            self.params.rates.selected_mutation_rate,
            self.params.rates.recombination_rate,
        )
        pdict["nregions"] = [fwdpy11.Region(0, 0.25, 1), fwdpy11.Region(0.5, 1, 1)]
        self.params = fwdpy11.ModelParams(**pdict)

    def test_neutral_mut_locations(self):
        fwdpy11.evolvets(self.rng, self.pop, self.params, 100)
        pos = [
            self.pop.tables.sites[i.site].position
            for i in self.pop.tables.mutations
            if i.neutral
        ]
        self.assertTrue(all([i < 0.25 or i >= 0.5 for i in pos]))


if __name__ == "__main__":
    unittest.main()
