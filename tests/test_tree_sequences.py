#
# Copyright (C) 2020 Kevin Thornton <krthornt@uci.edu>
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

import copy
import os
import pickle
import unittest
from collections import namedtuple

import fwdpy11
import fwdpy11.tskit_tools
import msprime
import numpy as np
import tskit


class Recorder(object):
    def __init__(self, seed, samplesize, timepoints):
        np.random.seed(seed)
        self.samplesize = samplesize
        self.timepoints = timepoints
        self.data = []

    def __call__(self, pop, recorder):
        if len(self.timepoints) > 0:
            if self.timepoints[0] == pop.generation:
                s = np.random.choice(pop.N, self.samplesize, replace=False)
                md = [copy.deepcopy(pop.diploid_metadata[i]) for i in s]
                self.data.append((pop.generation, md))
                recorder.assign(s)
                self.timepoints.pop(0)


def set_up_quant_trait_model(simlen=1.0):
    # TODO add neutral variants
    N = 1000
    demography = fwdpy11.DiscreteDemography()
    rho = 2.0
    # theta = 100.
    # nreps = 500
    # mu = theta/(4*N)
    r = rho / (4 * N)
    Opt = fwdpy11.Optimum
    GSSmo = fwdpy11.GSSmo(
        [Opt(when=0, optimum=0.0, VS=1.0), Opt(when=N, optimum=1.0, VS=1.0)]
    )
    a = fwdpy11.Additive(2.0, GSSmo)
    p = {
        "nregions": [],
        "sregions": [fwdpy11.GaussianS(0, 1, 1, 0.25)],
        "recregions": [fwdpy11.Region(0, 1, 1)],
        "rates": (0.0, 0.025, r),
        "gvalue": a,
        "prune_selected": False,
        "demography": demography,
        "simlen": np.rint(simlen * N).astype(int),
    }
    params = fwdpy11.ModelParams(**p)
    rng = fwdpy11.GSLrng(101 * 45 * 110 * 210)
    pop = fwdpy11.DiploidPopulation(N, 1.0)
    return params, rng, pop


def set_up_standard_pop_gen_model(simlen=1.0):
    """
    For this sort of model, when mutations fix, they are
    removed from the simulation, INCLUDING THE TREE
    SEQUENCES. The fact of their existence gets
    recorded in pop.fixations and pop.fixation_times
    """
    # TODO add neutral variants
    N = 1000
    rho = 1.0
    # theta = 100.
    # nreps = 500
    # mu = theta/(4*N)
    r = rho / (4 * N)

    a = fwdpy11.Multiplicative(2.0)
    pselected = 1e-3
    p = {
        "nregions": [],
        "sregions": [
            fwdpy11.GammaS(
                0, 1, 1.0 - pselected, mean=-5, shape_parameter=1, scaling=2 * N
            ),
            fwdpy11.ConstantS(0, 1, pselected, 1000, scaling=2 * N),
        ],
        "recregions": [fwdpy11.Region(0, 1, 1)],
        "rates": (0.0, 0.001, r),
        "gvalue": a,
        "prune_selected": True,
        "demography": fwdpy11.DiscreteDemography(),
        "simlen": np.rint(simlen * N).astype(int),
    }
    params = fwdpy11.ModelParams(**p)
    rng = fwdpy11.GSLrng(666 ** 2)
    pop = fwdpy11.DiploidPopulation(N, 1.0)
    return params, rng, pop


def mcounts_comparison_details(pop, counts, ts):
    for t in ts.trees():
        for m in t.mutations():
            pos = m.position
            # retrieve this mutation from the fwdpy11 table
            mr = None
            for i in pop.tables.mutations:
                if pop.tables.sites[i.site].position == pos:
                    mr = i
                    break
            assert pop.mutations[mr.key].pos == pos
            tskit_samples = [i for i in t.samples(m.node)]
            if len(tskit_samples) != counts[mr.key]:
                return False
    return True


def mcounts_comparison(pop, ts):
    if len(pop.tables.preserved_nodes) > 0:
        ts2 = ts.simplify(pop.tables.preserved_nodes)
        c = mcounts_comparison_details(pop, pop.mcounts_ancient_samples, ts2)
        if c is not True:
            return c

    ts2 = ts.simplify([i for i in range(2 * pop.N)])
    c = mcounts_comparison_details(pop, pop.mcounts, ts2)
    if c is not True:
        return c

    return True


def validate_mut_lookup_content(pop):
    ml = pop.mut_lookup
    ml_contents = []
    for i in ml:
        for j in ml[i]:
            ml_contents.append((i, j))
    ml_contents = sorted(ml_contents)

    mt_contents = []
    for mr in pop.tables.mutations:
        mt_contents.append((pop.tables.sites[mr.site].position, mr.key))
    mt_contents = sorted(mt_contents)
    return ml_contents == mt_contents


class TestBadGenomeLength(unittest.TestCase):
    def test_length_zero(self):
        with self.assertRaises(ValueError):
            fwdpy11.DiploidPopulation(100, 0.0)

    def test_negative_length(self):
        with self.assertRaises(ValueError):
            fwdpy11.DiploidPopulation(100, -1.0)


class TestTreeSequencesNoAncientSamplesKeepFixations(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.params, self.rng, self.pop = set_up_quant_trait_model(3.0)
        fwdpy11.evolvets(self.rng, self.pop, self.params, 100)
        assert (
            max(self.pop.mcounts) == 2 * self.pop.N
        ), "Nothing fixed, so test case is not helpful"

    def test_mut_lookup(self):
        self.assertEqual(len(self.pop.mut_lookup), len(self.pop.tables.mutations))
        self.assertEqual(len(self.pop.mut_lookup), len(self.pop.tables.sites))
        self.assertTrue(validate_mut_lookup_content(self.pop))

    def test_simplify(self):
        tables, idmap = fwdpy11.simplify(self.pop, [i for i in range(10)])
        for i in range(10):
            self.assertTrue(idmap[i] != fwdpy11.NULL_NODE)

    def test_simplify_tables(self):
        tables, idmap = fwdpy11.simplify_tables(self.pop.tables, [i for i in range(10)])
        for i in range(10):
            self.assertTrue(idmap[i] != fwdpy11.NULL_NODE)

    def test_simplify_numpy_array(self):
        tables, idmap = fwdpy11.simplify(self.pop, np.array([i for i in range(10)]))
        for i in range(10):
            self.assertTrue(idmap[i] != fwdpy11.NULL_NODE)

    def test_simplify_tables_numpy_array(self):
        tables, idmap = fwdpy11.simplify_tables(
            self.pop.tables, np.array([i for i in range(10)])
        )
        for i in range(10):
            self.assertTrue(idmap[i] != fwdpy11.NULL_NODE)

    def test_dump_to_tskit(self):
        import tskit

        dumped_ts = self.pop.dump_tables_to_tskit()
        self.assertEqual(len(dumped_ts.tables.nodes), len(self.pop.tables.nodes))
        self.assertEqual(len(dumped_ts.tables.edges), len(self.pop.tables.edges))
        self.assertEqual(
            len(dumped_ts.tables.mutations), len(self.pop.tables.mutations)
        )
        eview = np.array(self.pop.tables.edges, copy=False)
        self.assertEqual(eview["parent"].sum(), dumped_ts.tables.edges.parent.sum())
        self.assertEqual(eview["child"].sum(), dumped_ts.tables.edges.child.sum())
        self.assertEqual(eview["left"].sum(), dumped_ts.tables.edges.left.sum())
        self.assertEqual(eview["right"].sum(), dumped_ts.tables.edges.right.sum())
        tv = fwdpy11.TreeIterator(self.pop.tables, [i for i in range(2 * self.pop.N)])
        tt_fwd = 0
        for t in tv:
            tt_fwd += t.total_time(self.pop.tables.nodes)
        tt_tskit = 0
        for t in dumped_ts.trees():
            tt_tskit += t.get_total_branch_length()
        self.assertEqual(tt_fwd, tt_tskit)

        for i, j in enumerate(self.pop.diploid_metadata):
            self.assertEqual(dumped_ts.tables.nodes.individual[j.nodes[0]], i)
            self.assertEqual(dumped_ts.tables.nodes.individual[j.nodes[1]], i)

        # Check that our individual flags, added in 0.8.0, are okay
        tskit_alive_nodes = (
            dumped_ts.tables.individuals.flags & fwdpy11.tskit_tools.INDIVIDUAL_IS_ALIVE
        )
        self.assertTrue(len((tskit_alive_nodes > 0)), len(self.pop.diploid_metadata))
        for i, j in enumerate(self.pop.diploid_metadata):
            self.assertTrue(
                (
                    dumped_ts.tables.individuals.flags[i]
                    & fwdpy11.tskit_tools.INDIVIDUAL_IS_PRESERVED
                )
                == 0
            )
            self.assertTrue(
                (
                    dumped_ts.tables.individuals.flags[i]
                    & fwdpy11.tskit_tools.INDIVIDUAL_IS_FIRST_GENERATION
                )
                == 0
            )
            self.assertTrue(
                (
                    dumped_ts.tables.individuals.flags[i]
                    & fwdpy11.tskit_tools.INDIVIDUAL_IS_ALIVE
                )
                > 0
            )
            for n in j.nodes:
                self.assertTrue(
                    (dumped_ts.tables.nodes.flags[n] & tskit.NODE_IS_SAMPLE) > 0
                )

        # Now, we make sure that the metadata can
        # be decoded
        md = []
        for i, j in enumerate(self.pop.diploid_metadata):
            d = dumped_ts.tables.individuals[i].metadata
            md.append(d)
            self.assertEqual(j.g, d["g"])
            self.assertEqual(j.w, d["w"])
            self.assertEqual(j.e, d["e"])
            self.assertEqual(j.label, d["label"])
            self.assertEqual(j.sex, d["sex"])
            self.assertEqual(j.deme, d["deme"])
            for k, l in zip(j.parents, d["parents"]):
                self.assertEqual(k, l)
            for k, l in zip(j.geography, d["geography"]):
                self.assertEqual(k, l)

        # Test that we can go backwards from node table to individuals
        samples = np.where(dumped_ts.tables.nodes.flags == tskit.NODE_IS_SAMPLE)[0]
        self.assertEqual(len(samples), 2 * self.pop.N)
        for i in samples[::2]:
            ind = i // 2
            d = md[ind]
            fwdpy11_md = self.pop.diploid_metadata[ind]
            self.assertEqual(fwdpy11_md.g, d["g"])
            self.assertEqual(fwdpy11_md.w, d["w"])
            self.assertEqual(fwdpy11_md.e, d["e"])
            self.assertEqual(fwdpy11_md.label, d["label"])
            self.assertEqual(fwdpy11_md.sex, d["sex"])
            self.assertEqual(fwdpy11_md.deme, d["deme"])
            for k, l in zip(fwdpy11_md.parents, d["parents"]):
                self.assertEqual(k, l)
            for k, l in zip(fwdpy11_md.geography, d["geography"]):
                self.assertEqual(k, l)

        for idx in range(len(self.pop.tables.mutations)):
            i = self.pop.tables.mutations[idx]
            j = dumped_ts.tables.mutations[idx].site
            d = dumped_ts.tables.mutations[idx].metadata
            self.assertEqual(i.key, d["key"])
            site = dumped_ts.tables.sites[j]
            m = self.pop.mutations[d["key"]]
            self.assertEqual(site.position, m.pos)
            self.assertEqual(d["s"], m.s)
            self.assertEqual(d["h"], m.h)
            self.assertEqual(d["label"], m.label)
            self.assertEqual(d["neutral"], m.neutral)

        self.assertEqual(mcounts_comparison(self.pop, dumped_ts), True)

    def test_dump_to_tskit_parameters_round_trip(self):
        import json
        import tskit

        dumped_ts = self.pop.dump_tables_to_tskit(
            parameters={"model": str(self.params)}
        )
        p = json.loads(dumped_ts.provenance(0).record)
        try:
            tskit.validate_provenance(p)
        except:  # NOQA
            self.fail("unexpected exception")

        array = np.array  # NOQA
        mp = eval(p["parameters"]["model"])
        self.assertTrue(mp == self.params)

    def test_TreeIterator(self):
        # The first test ensures that TreeIterator
        # simply holds a reference to the input tables,
        # rather than a (deep) copy, which would have
        # a different address
        tv = fwdpy11.TreeIterator(
            self.pop.tables, [i for i in range(2 * self.pop.N)], True, 0, 1
        )
        self.assertTrue(tv.tables is self.pop.tables)

        with self.assertRaises(ValueError):
            tv = fwdpy11.TreeIterator(
                self.pop.tables, [i for i in range(2 * self.pop.N)], True, 1, 0
            )
        with self.assertRaises(ValueError):
            tv = fwdpy11.TreeIterator(
                self.pop.tables, [i for i in range(2 * self.pop.N)], False, 1, 0
            )

        for i in np.arange(0.0, 1.0, 0.1):
            tv = fwdpy11.TreeIterator(
                self.pop.tables, [i for i in range(2 * self.pop.N)], True, i, i + 0.1
            )
            for ti in tv:
                a = ti.left < i + 0.1
                b = i < ti.right
                self.assertTrue(a and b)

            tv = fwdpy11.TreeIterator(
                self.pop.tables, [i for i in range(2 * self.pop.N)], False, i, i + 0.1
            )
            for ti in tv:
                a = ti.left < i + 0.1
                b = i < ti.right
                self.assertTrue(a and b)

    def test_TreeIterator_iterate_sites(self):
        # TODO: need test of empty tree sequence
        # and tree sequence where mutations aren't
        # on every tree
        tv = fwdpy11.TreeIterator(self.pop.tables, [i for i in range(2 * self.pop.N)])
        nsites_visited = 0
        for tree in tv:
            for s in tree.sites():
                self.assertTrue(s.position >= tree.left)
                self.assertTrue(s.position < tree.right)
                nsites_visited += 1
        self.assertEqual(len(self.pop.tables.sites), nsites_visited)

        site_table = np.array(self.pop.tables.sites, copy=False)
        for i in np.arange(0.0, 1.0, 0.1):
            tv = fwdpy11.TreeIterator(
                self.pop.tables, [i for i in range(2 * self.pop.N)], False, i, i + 0.1
            )
            nsites_visited = 0
            idx = np.where(
                (site_table["position"] >= i) & (site_table["position"] < i + 0.1)
            )[0]
            nsites_in_interval = len(idx)
            for tree in tv:
                for s in tree.sites():
                    self.assertTrue(s.position >= tree.left)
                    self.assertTrue(s.position < tree.right)
                    nsites_visited += 1
            self.assertEqual(nsites_visited, nsites_in_interval)

    def test_TreeIterator_iterate_mutations(self):
        tv = fwdpy11.TreeIterator(self.pop.tables, [i for i in range(2 * self.pop.N)])
        sites = np.array(self.pop.tables.sites, copy=False)
        nsites_visited = 0
        for tree in tv:
            for m in tree.mutations():
                self.assertTrue(sites["position"][m.site] >= tree.left)
                self.assertTrue(sites["position"][m.site] < tree.right)
                nsites_visited += 1
        self.assertEqual(len(self.pop.tables.sites), nsites_visited)

        for i in np.arange(0.0, 1.0, 0.1):
            tv = fwdpy11.TreeIterator(
                self.pop.tables, [i for i in range(2 * self.pop.N)], False, i, i + 0.1
            )
            nsites_visited = 0
            idx = np.where((sites["position"] >= i) & (sites["position"] < i + 0.1))[0]
            nsites_in_interval = len(idx)
            for tree in tv:
                for m in tree.mutations():
                    self.assertTrue(sites["position"][m.site] >= tree.left)
                    self.assertTrue(sites["position"][m.site] < tree.right)
                    nsites_visited += 1
            self.assertEqual(nsites_visited, nsites_in_interval)

    def test_leaf_counts_vs_mcounts(self):
        tv = fwdpy11.TreeIterator(self.pop.tables, [i for i in range(2 * self.pop.N)])
        mv = np.array(self.pop.tables.mutations, copy=False)
        muts = self.pop.mutations_ndarray
        p = muts["pos"]
        for t in tv:
            l, r = t.left, t.right
            mt = [i for i in mv if p[i[1]] >= l and p[i[1]] < r]
            for i in mt:
                self.assertEqual(t.leaf_counts(i[0]), self.pop.mcounts[i[1]])

    def test_simplify_to_sample(self):
        """
        Simplify to a sample using fwdpy11 and
        tskit, then test that total time on output
        is the same from both sources and that
        the mutation tables contain the same
        positions after simplification.
        """
        dumped_ts = self.pop.dump_tables_to_tskit()
        tt = 0.0
        for i in self.pop.tables.nodes:
            tt += i.time
        samples = np.arange(0, 2 * self.pop.N, 50, dtype=np.int32)
        mspts = dumped_ts.simplify(samples=samples.tolist())
        fp11ts, idmap = fwdpy11.simplify(self.pop, samples)
        for i in range(len(fp11ts.edges)):
            self.assertTrue(fp11ts.edges[i].parent < len(fp11ts.nodes))
            self.assertTrue(fp11ts.edges[i].child < len(fp11ts.nodes))
        for s in samples:
            self.assertEqual(fp11ts.nodes[idmap[s]].time, self.pop.generation)
        tt_fwd = 0.0
        tv = fwdpy11.TreeIterator(fp11ts, [i for i in range(len(samples))])
        for t in tv:
            tt_fwd += t.total_time(fp11ts.nodes)
        tt_tskit = 0.0
        for t in mspts.trees():
            tt_tskit += t.get_total_branch_length()
        self.assertEqual(tt_fwd, tt_tskit)

        self.assertEqual(len(fp11ts.mutations), len(mspts.tables.mutations))
        fp11_pos = np.array([self.pop.mutations[i.key].pos for i in fp11ts.mutations])
        fp11_pos = np.sort(fp11_pos)
        msp_pos = np.sort(mspts.tables.sites.position)
        self.assertTrue(np.array_equal(fp11_pos, msp_pos))

    def test_genotype_matrix(self):
        """
        Make data matrix objects from the tree sequences
        and compare their contents to those of tskit
        as well as to an explicit calculation of mutation counts.
        """
        dm = fwdpy11.data_matrix_from_tables(
            self.pop.tables, [i for i in range(2 * self.pop.N)], False, True, True
        )
        sa = np.array(dm.selected)
        cs = np.sum(sa, axis=1)
        dumped_ts = self.pop.dump_tables_to_tskit()
        mm = dumped_ts.genotype_matrix()
        mc = np.sum(mm, axis=1)
        ec = np.zeros(len(self.pop.mutations), dtype=np.uint32)
        for d in self.pop.diploids:
            for m in self.pop.haploid_genomes[d.first].smutations:
                ec[m] += 1
            for m in self.pop.haploid_genomes[d.second].smutations:
                ec[m] += 1
        for i, j, k in zip(self.pop.tables.mutations, cs, mc):
            self.assertEqual(self.pop.mcounts[i.key], ec[i.key])
            self.assertEqual(ec[i.key], j)
            self.assertEqual(ec[i.key], k)
        self.assertTrue(np.array_equal(ec, np.array(self.pop.mcounts)))
        self.assertEqual(ec.sum(), cs.sum())
        self.assertEqual(ec.sum(), mc.sum())
        self.assertTrue(np.array_equal(sa, mm))
        self.assertTrue(np.array_equal(cs, mc))

    def test_genotype_matrix_ranges(self):
        dm = fwdpy11.data_matrix_from_tables(
            self.pop.tables, [i for i in range(2 * self.pop.N)], False, True
        )
        spos = np.array(dm.selected.positions)
        for i in np.arange(0, self.pop.tables.genome_length, 0.1):
            dmi = fwdpy11.data_matrix_from_tables(
                self.pop.tables,
                [i for i in range(2 * self.pop.N)],
                record_neutral=False,
                record_selected=True,
                begin=i,
                end=i + 0.1,
            )
            w = np.where((spos >= i) & (spos < i + 0.1))[0]
            self.assertTrue(np.array_equal(spos[w], np.array(dmi.selected.positions)))

    def test_VariantIterator(self):
        """
        Test VariantIterator by asserting
        that sum of genotypes equal values in
        the corresponding DataMatrix and
        those in pop.mcounts
        """
        dm = fwdpy11.data_matrix_from_tables(
            self.pop.tables, [i for i in range(2 * self.pop.N)], False, True, True
        )
        sa = np.array(dm.selected)
        cs = np.sum(sa, axis=1)
        i = 0
        vi = fwdpy11.VariantIterator(
            self.pop.tables, [i for i in range(2 * self.pop.N)]
        )
        for v in vi:
            c = self.pop.mcounts[self.pop.tables.mutations[i].key]
            self.assertEqual(c, cs[i])
            self.assertEqual(c, v.genotypes.sum())
            i += 1
        mc = np.array(self.pop.mcounts)
        self.assertEqual(i, len(np.where(mc > 0)[0]))
        self.assertEqual(i, len(self.pop.tables.mutations))

    def test_VariantIteratorFromPopulation(self):
        dm = fwdpy11.data_matrix_from_tables(
            self.pop.tables, [i for i in range(2 * self.pop.N)], False, True, True
        )
        sa = np.array(dm.selected)
        cs = np.sum(sa, axis=1)
        i = 0
        vi = fwdpy11.VariantIterator(self.pop)
        for v in vi:
            c = self.pop.mcounts[self.pop.tables.mutations[i].key]
            self.assertEqual(c, cs[i])
            self.assertEqual(c, v.genotypes.sum())
            i += 1
        mc = np.array(self.pop.mcounts)
        self.assertEqual(i, len(np.where(mc > 0)[0]))
        self.assertEqual(i, len(self.pop.tables.mutations))

    def test_VariantIteratorBeginEnd(self):
        for i in np.arange(0, self.pop.tables.genome_length, 0.1):
            vi = fwdpy11.VariantIterator(
                self.pop.tables, [i for i in range(2 * self.pop.N)], i, i + 0.1
            )
            nm = len(
                [
                    j
                    for j in self.pop.tables.mutations
                    if self.pop.mutations[j.key].pos >= i
                    and self.pop.mutations[j.key].pos < i + 0.1
                ]
            )
            nseen = 0
            for v in vi:
                r = v.records[0]
                self.assertTrue(self.pop.mutations[r.key].pos >= i)
                self.assertTrue(self.pop.mutations[r.key].pos < i + 0.1)
                nseen += 1
            self.assertEqual(nm, nseen)

        # test bad start/stop
        with self.assertRaises(ValueError):
            vi = fwdpy11.VariantIterator(
                self.pop.tables, [i for i in range(2 * self.pop.N)], begin=0.5, end=0.25
            )

    def test_count_mutations(self):
        mc = fwdpy11.count_mutations(self.pop, [i for i in range(2 * self.pop.N)])
        pmc = np.array(self.pop.mcounts)
        self.assertTrue(np.array_equal(mc, pmc))

    def test_binary_round_trip(self):
        ofile = "poptest_no_ancient_preserve_fixations.bin"
        self.pop.dump_to_file(ofile)
        pop2 = fwdpy11.DiploidPopulation.load_from_file(ofile)
        self.assertTrue(self.pop == pop2)
        if os.path.exists(ofile):
            os.remove(ofile)

    def test_fast_pickling(self):
        p = pickle.dumps(self.pop, -1)
        up = pickle.loads(p)
        self.assertTrue(self.pop == up)

    def test_slow_pickling_to_file(self):
        ofile = "poptest_no_ancient_preserve_fixations.pickle"
        with open(ofile, "wb") as f:
            self.pop.pickle_to_file(f)
        with open(ofile, "rb") as f:
            pop2 = fwdpy11.DiploidPopulation.load_from_pickle_file(f)
        self.assertTrue(self.pop == pop2)
        if os.path.exists(ofile):
            os.remove(ofile)


class TestTreeSequencesWithAncientSamplesKeepFixations(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.params, self.rng, self.pop = set_up_quant_trait_model(3.0)
        self.stimes = [i for i in range(1, 101)]
        self.recorder = Recorder(42, 10, self.stimes)
        fwdpy11.evolvets(self.rng, self.pop, self.params, 100, self.recorder)
        assert (
            max(self.pop.mcounts) == 2 * self.pop.N
        ), "Nothing fixed, so test case is not helpful"

    def test_mut_lookup(self):
        self.assertEqual(len(self.pop.mut_lookup), len(self.pop.tables.mutations))
        self.assertEqual(len(self.pop.mut_lookup), len(self.pop.tables.sites))
        self.assertTrue(validate_mut_lookup_content(self.pop))

    def test_serialization_via_file_round_trip(self):
        with open("test_ts_serialization_via_file.pickle", "wb") as f:
            self.pop.pickle_to_file(f)

        with open("test_ts_serialization_via_file.pickle", "rb") as f:
            up = fwdpy11.DiploidPopulation.load_from_pickle_file(f)
        self.assertTrue(self.pop == up)
        self.assertTrue(self.pop.mut_lookup == up.mut_lookup)

    def test_pickling_round_trip(self):
        p = pickle.dumps(self.pop, -1)
        up = pickle.loads(p)
        self.assertTrue(self.pop == up)
        self.assertTrue(self.pop.mut_lookup == up.mut_lookup)

    def test_count_mutations_preserved_samples(self):
        mc = fwdpy11.count_mutations(self.pop, self.pop.tables.preserved_nodes)
        pmc = np.array(self.pop.mcounts_ancient_samples)
        self.assertTrue(np.array_equal(mc, pmc))

    def test_ancient_sample_times(self):
        times = []
        for t, n, md in self.pop.sample_timepoints(False):
            times.append(int(t))
        self.assertEqual(times, [i for i in range(1, 101)])

    def test_VariantIteratorFromPreservedSamples(self):
        n = np.array(self.pop.tables.nodes)
        pn = np.array(self.pop.tables.preserved_nodes)
        at = n["time"][pn]
        for u in np.unique(at):
            n = pn[np.where(at == u)[0]]
            vi = fwdpy11.VariantIterator(self.pop.tables, n)
            for variant in vi:
                k = variant.records[0]
                self.assertNotEqual(k.node, fwdpy11.NULL_NODE)
                self.assertNotEqual(k.key, np.iinfo(np.uint64).max)

    def test_mcounts_from_ancient_samples(self):
        vi = fwdpy11.VariantIterator(self.pop.tables, self.pop.tables.preserved_nodes)
        for v in vi:
            k = v.records[0]
            self.assertEqual(self.pop.mcounts_ancient_samples[k.key], v.genotypes.sum())

    def test_sample_traverser(self):
        timepoints = [i for i in range(1, 101)]
        amd = np.array(self.pop.ancient_sample_metadata, copy=False)
        n = np.array(self.pop.tables.nodes)
        at = n["time"][amd["nodes"][:, 0]]
        pn = np.array(self.pop.tables.preserved_nodes)
        self.assertEqual(2 * len(amd), len(pn))
        # j contains (time, nodes, metadata).  The metadata
        # also contain nodes
        for i, j, d in zip(
            timepoints, self.pop.sample_timepoints(False), self.recorder.data
        ):
            self.assertEqual(i, d[0])
            self.assertEqual(i, j[0])
            self.assertTrue(all(n["time"][j[1]] == i))
            mdn = j[2]["nodes"].flatten()
            self.assertTrue(np.array_equal(j[1], mdn))
            idx = np.where(at == i)[0]
            self.assertTrue(np.array_equal(j[2], amd[idx]))
            # Compare metadata stored in recorder
            # to that preserved internally by the simulation
            for k, l in zip(idx, d[1]):
                self.assertEqual(amd["g"][k], l.g)
                self.assertEqual(amd["w"][k], l.w)
                self.assertEqual(amd["e"][k], l.e)
                self.assertEqual(amd["label"][k], l.label)
                self.assertEqual(amd["parents"][k][0], l.parents[0])
                self.assertEqual(amd["parents"][k][1], l.parents[1])

            # Extract out the nodes from preserved_nodes
            idx = np.where(n["time"][self.pop.tables.preserved_nodes] == i)[0]
            self.assertTrue(
                np.array_equal(np.array(self.pop.tables.preserved_nodes)[idx], j[1])
            )

    def test_binary_round_trip(self):
        ofile = "poptest_with_ancient_preserve_fixations.bin"
        self.pop.dump_to_file(ofile)
        pop2 = fwdpy11.DiploidPopulation.load_from_file(ofile)
        self.assertTrue(self.pop == pop2)
        self.assertTrue(self.pop.mut_lookup == pop2.mut_lookup)
        if os.path.exists(ofile):
            os.remove(ofile)

    def test_fast_pickling(self):
        p = pickle.dumps(self.pop, -1)
        up = pickle.loads(p)
        self.assertTrue(self.pop == up)

    def test_slow_pickling_to_file(self):
        ofile = "poptest_with_ancient_preserve_fixations.pickle"
        with open(ofile, "wb") as f:
            self.pop.pickle_to_file(f)
        with open(ofile, "rb") as f:
            pop2 = fwdpy11.DiploidPopulation.load_from_pickle_file(f)
        self.assertTrue(self.pop == pop2)
        if os.path.exists(ofile):
            os.remove(ofile)

    def test_dump_to_tskit(self):
        dumped_ts = self.pop.dump_tables_to_tskit()
        self.assertEqual(mcounts_comparison(self.pop, dumped_ts), True)

        metadata = [i for i in self.pop.diploid_metadata]
        metadata.extend([i for i in self.pop.ancient_sample_metadata])

        self.assertTrue(len(metadata) > len(self.pop.diploid_metadata))
        for i, j in enumerate(metadata):
            self.assertEqual(dumped_ts.tables.nodes.individual[j.nodes[0]], i)
            self.assertEqual(dumped_ts.tables.nodes.individual[j.nodes[1]], i)

        # Check that our individual flags, added in 0.8.0, are okay
        tskit_alive_nodes = (
            dumped_ts.tables.individuals.flags & fwdpy11.tskit_tools.INDIVIDUAL_IS_ALIVE
        )
        self.assertTrue(len((tskit_alive_nodes > 0)), len(self.pop.diploid_metadata))
        for i, j in enumerate(self.pop.ancient_sample_metadata):
            for n in j.nodes:
                self.assertTrue(
                    (dumped_ts.tables.nodes.flags[n] & tskit.NODE_IS_SAMPLE) > 0
                )
            self.assertTrue(
                (
                    dumped_ts.tables.individuals.flags[
                        i + len(self.pop.diploid_metadata)
                    ]
                    & fwdpy11.tskit_tools.INDIVIDUAL_IS_PRESERVED
                )
                > 1
            )
            self.assertTrue(
                (
                    dumped_ts.tables.individuals.flags[
                        i + len(self.pop.diploid_metadata)
                    ]
                    & fwdpy11.tskit_tools.INDIVIDUAL_IS_FIRST_GENERATION
                )
                == 0
            )
            self.assertTrue(
                (
                    dumped_ts.tables.individuals.flags[
                        i + len(self.pop.diploid_metadata)
                    ]
                    & fwdpy11.tskit_tools.INDIVIDUAL_IS_ALIVE
                )
                == 0
            )


class TestMutationCounts(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.params, self.rng, self.pop = set_up_quant_trait_model()
        pdict = self.params.asdict()
        pdict["simlen"] = 200
        self.params = fwdpy11.ModelParams(**pdict)
        self.pop2 = copy.deepcopy(self.pop)
        self.rng2 = fwdpy11.GSLrng(101 * 45 * 110 * 210)
        params2 = copy.deepcopy(self.params)
        fwdpy11.evolvets(self.rng, self.pop, self.params, 100)
        fwdpy11.evolvets(
            self.rng2, self.pop2, params2, 100, remove_extinct_variants=False
        )

    def test_mut_lookup(self):
        self.assertEqual(len(self.pop.mut_lookup), len(self.pop.tables.mutations))
        self.assertEqual(len(self.pop.mut_lookup), len(self.pop.tables.sites))
        self.assertTrue(validate_mut_lookup_content(self.pop))

    def test_mutation_counts(self):
        self.assertEqual(sum(self.pop.mcounts), sum(self.pop2.mcounts))

    def test_mutation_counts_details(self):
        mc = [i for i in self.pop2.mcounts if i > 0]
        self.assertTrue(all([i == j for i, j in zip(self.pop.mcounts, mc)]) is True)

    def test_mutation_counts_from_genomes(self):
        """
        Removing extinct genomes happens after
        the final mutation count, so let's
        make sure we can reconstruct the counts
        from the remapped genomes
        """
        mc = [0] * len(self.pop.mutations)
        mc2 = [0] * len(self.pop2.mutations)
        for g in self.pop.haploid_genomes:
            if g.n > 0:
                for k in g.smutations:
                    mc[k] += g.n
        for g in self.pop2.haploid_genomes:
            if g.n > 0:
                for k in g.smutations:
                    mc2[k] += g.n
        self.assertTrue(all([i == j for i, j in zip(mc, self.pop.mcounts)]) is True)
        self.assertTrue(all([i == j for i, j in zip(mc2, self.pop2.mcounts)]) is True)


class TestUpdateTiming(unittest.TestCase):
    """
    Motivated by issue 437.
    This test case constructs a situation
    where we know exactly when a change should
    happen and what its outcome should be.
    If the optimum shifts at generation X,
    then the current parents should experience
    reduced fitness at that time, which we can
    measure by their mean fitness, which is trivial
    if there are no mutations.
    """

    def testTiming(self):
        """
        Make sure that timings occur when we expect.
        """
        from test_tree_sequences import set_up_quant_trait_model

        params, rng, pop = set_up_quant_trait_model(1.1)
        pdict = params.asdict()
        pdict["rates"] = (0, 0, 0)
        params = fwdpy11.ModelParams(**pdict)

        class Recorder(object):
            def __init__(self):
                self.data = []

            def __call__(self, pop, sampler):
                if pop.generation >= 990:
                    md = np.array(pop.diploid_metadata, copy=False)
                    self.data.append((pop.generation, md["w"].mean()))

        r = Recorder()
        for i in r.data:
            if i[0] < pop.N:
                self.assertAlmostEqual(i[1], 1.0)
            else:
                self.assertAlmostEqual(i[1], np.exp(-0.50))


class TestTreeSequencesNoAncientSamplesPruneFixations(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.params, self.rng, self.pop = set_up_standard_pop_gen_model()
        fwdpy11.evolvets(
            self.rng, self.pop, self.params, 100, track_mutation_counts=True
        )
        assert len(self.pop.fixations) > 0, "Nothing fixed, so test case is not helpful"

    def test_mut_lookup(self):
        self.assertEqual(len(self.pop.mut_lookup), len(self.pop.tables.mutations))
        self.assertEqual(len(self.pop.mut_lookup), len(self.pop.tables.sites))
        self.assertTrue(validate_mut_lookup_content(self.pop))

    def test_max_mcounts(self):
        self.assertTrue(max(self.pop.mcounts) < 2 * self.pop.N)

    def test_mutation_table_contents(self):
        self.assertEqual(len(self.pop.mcounts), len(self.pop.mcounts_ancient_samples))
        for m in self.pop.tables.mutations:
            self.assertTrue(m.key < len(self.pop.mutations))
            self.assertTrue(self.pop.mcounts[m.key] < 2 * self.pop.N)

        for f in self.pop.fixations:
            m = None
            for i in self.pop.tables.mutations:
                if self.pop.tables.sites[i.site].position == f.pos:
                    m = i
                    break
            if m is not None:
                # Then the mutation position got re-used, and
                # so the mutations must have different origin times
                # NOTE: this probably has not happened!
                self.assertNotEqual(f.g, self.pop.mutations[m.key].g)

    def test_dump_to_tskit(self):
        ts = self.pop.dump_tables_to_tskit()
        self.assertEqual(mcounts_comparison(self.pop, ts), True)

    def test_binary_round_trip(self):
        ofile = "poptest_no_ancient_prune_fixations.bin"
        self.pop.dump_to_file(ofile)
        pop2 = fwdpy11.DiploidPopulation.load_from_file(ofile)
        self.assertTrue(self.pop == pop2)
        if os.path.exists(ofile):
            os.remove(ofile)

    def test_fast_pickling(self):
        p = pickle.dumps(self.pop, -1)
        up = pickle.loads(p)
        self.assertTrue(self.pop == up)

    def test_slow_pickling_to_file(self):
        ofile = "poptest_no_ancient_prune_fixations.pickle"
        with open(ofile, "wb") as f:
            self.pop.pickle_to_file(f)
        with open(ofile, "rb") as f:
            pop2 = fwdpy11.DiploidPopulation.load_from_pickle_file(f)
        self.assertTrue(self.pop == pop2)
        if os.path.exists(ofile):
            os.remove(ofile)

    def test_genotype_matrix(self):
        dm = fwdpy11.data_matrix_from_tables(
            self.pop.tables, [i for i in range(2 * self.pop.N)], False, True, True
        )
        rc = np.sum(dm.selected, axis=1)
        index = [i for i in range(len(self.pop.mcounts))]
        index = sorted(index, key=lambda x: self.pop.mutations[x].pos)
        mc = [self.pop.mcounts[i] for i in index if self.pop.mcounts[i] > 0]
        self.assertTrue(np.array_equal(rc, np.array(mc)))


class TestTreeSequencesWithAncientSamplesPruneFixations(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.params, self.rng, self.pop = set_up_standard_pop_gen_model()
        self.stimes = [i for i in range(1, 101)]
        self.recorder = Recorder(42, 10, self.stimes)
        fwdpy11.evolvets(
            self.rng,
            self.pop,
            self.params,
            100,
            self.recorder,
            track_mutation_counts=True,
        )
        assert len(self.pop.fixations) > 0, "Nothing fixed, so test case is not helpful"

    def test_mut_lookup(self):
        self.assertEqual(len(self.pop.mut_lookup), len(self.pop.tables.mutations))
        self.assertEqual(len(self.pop.mut_lookup), len(self.pop.tables.sites))
        self.assertTrue(validate_mut_lookup_content(self.pop))

    def test_mutation_table_contents(self):
        self.assertEqual(len(self.pop.mcounts), len(self.pop.mcounts_ancient_samples))
        for m in self.pop.tables.mutations:
            self.assertTrue(m.key < len(self.pop.mutations))
            if self.pop.mcounts[m.key] < 2 * self.pop.N:
                is_fixed = False
            else:
                is_fixed = True
            if self.pop.mcounts_ancient_samples[m.key] != 0:
                is_found_ancient = True
            else:
                is_found_ancient = False
            self.assertTrue(
                is_fixed is False or (is_fixed is True and is_found_ancient is True)
            )

        for f in self.pop.fixations:
            m = None
            for i in self.pop.tables.mutations:
                if self.pop.tables.sites[i.site].position == f.pos:
                    m = i
                    break
            if m is not None:
                # Then the mutation position got re-used, and
                # so the mutations must have different origin times
                # NOTE: this probably has not happened!
                self.assertNotEqual(f.g, self.pop.mutations[m.key].g)

    def test_dump_to_tskit(self):
        ts = self.pop.dump_tables_to_tskit()
        self.assertEqual(mcounts_comparison(self.pop, ts), True)

    def test_binary_round_trip(self):
        ofile = "poptest_with_ancient_prune_fixations.bin"
        self.pop.dump_to_file(ofile)
        pop2 = fwdpy11.DiploidPopulation.load_from_file(ofile)
        self.assertTrue(self.pop == pop2)
        if os.path.exists(ofile):
            os.remove(ofile)

    def test_fast_pickling(self):
        p = pickle.dumps(self.pop, -1)
        up = pickle.loads(p)
        self.assertTrue(self.pop == up)

    def test_slow_pickling_to_file(self):
        ofile = "poptest_with_ancient_prune_fixations.pickle"
        with open(ofile, "wb") as f:
            self.pop.pickle_to_file(f)
        with open(ofile, "rb") as f:
            pop2 = fwdpy11.DiploidPopulation.load_from_pickle_file(f)
        self.assertTrue(self.pop == pop2)
        if os.path.exists(ofile):
            os.remove(ofile)

    def test_genotype_matrix(self):
        dm = fwdpy11.data_matrix_from_tables(
            self.pop.tables, [i for i in range(2 * self.pop.N)], False, True, True
        )
        rc = np.sum(dm.selected, axis=1)
        index = [i for i in range(len(self.pop.mcounts))]
        index = sorted(index, key=lambda x: self.pop.mutations[x].pos)
        mc = [self.pop.mcounts[i] for i in index if self.pop.mcounts[i] > 0]
        self.assertTrue(np.array_equal(rc, np.array(mc)))


class TestSimplificationInterval(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.N = 1000
        self.rho = 1.0
        self.theta = 100.0
        self.nreps = 500
        self.mu = self.theta / (4 * self.N)
        self.r = self.rho / (4 * self.N)
        self.GSS = fwdpy11.GSS(VS=1.0, optimum=0.0)
        a = fwdpy11.Additive(2.0, self.GSS)
        self.p = {
            "nregions": [],
            "sregions": [fwdpy11.GaussianS(0, 1, 1, 0.25)],
            "recregions": [fwdpy11.Region(0, 1, 1)],
            "rates": (0.0, 0.025, self.r),
            "gvalue": a,
            "prune_selected": False,
            "demography": fwdpy11.DiscreteDemography(),
            "simlen": 100,
        }
        self.params = fwdpy11.ModelParams(**self.p)
        self.rng = fwdpy11.GSLrng(101 * 45 * 110 * 210)
        self.pop = fwdpy11.DiploidPopulation(self.N, 1.0)
        self.recorder = fwdpy11.RandomAncientSamples(
            seed=42, samplesize=10, timepoints=[i for i in range(1, 101)]
        )

    def testEvolve(self):
        # TODO: actually test something here :)
        fwdpy11.evolvets(self.rng, self.pop, self.params, 1, self.recorder)
        samples = [i for i in range(2 * self.pop.N)] + self.pop.tables.preserved_nodes
        vi = fwdpy11.TreeIterator(self.pop.tables, samples)

    def test_mut_lookup(self):
        self.assertEqual(len(self.pop.mut_lookup), len(self.pop.tables.mutations))
        self.assertEqual(len(self.pop.mut_lookup), len(self.pop.tables.sites))
        self.assertTrue(validate_mut_lookup_content(self.pop))


class TestFixationPreservation(unittest.TestCase):
    def testQtraitSim(self):
        N = 1000
        rho = 1.0
        r = rho / (4 * N)

        GSS = fwdpy11.GSS(VS=1.0, optimum=1.0)
        a = fwdpy11.Additive(2.0, GSS)
        p = {
            "nregions": [],
            "sregions": [fwdpy11.GaussianS(0, 1, 1, 0.25)],
            "recregions": [fwdpy11.Region(0, 1, 1)],
            "rates": (0.0, 0.005, r),
            "gvalue": a,
            "prune_selected": False,
            "demography": fwdpy11.DiscreteDemography(),
            "simlen": 3 * N,
        }
        params = fwdpy11.ModelParams(**p)
        rng = fwdpy11.GSLrng(101 * 45 * 110 * 210)
        pop = fwdpy11.DiploidPopulation(N, 1.0)
        fwdpy11.evolvets(rng, pop, params, 100, track_mutation_counts=True)
        mc = fwdpy11.count_mutations(
            pop.tables, pop.mutations, [i for i in range(2 * pop.N)]
        )
        assert len(pop.fixations) > 0, "Test is meaningless without fixations"
        fixations = np.where(mc == 2 * pop.N)[0]
        self.assertEqual(len(fixations), len(pop.fixations))

        # Brute-force calculation of fixations
        brute_force = np.zeros(len(pop.mutations), dtype=np.int32)
        for g in pop.haploid_genomes:
            if g.n > 0:
                for k in g.smutations:
                    brute_force[k] += g.n

        self.assertTrue(np.array_equal(brute_force, mc))
        self.assertTrue(np.array_equal(brute_force, pop.mcounts))

    def testPopGenSim(self):
        N = 1000
        rho = 1.0
        r = rho / (4 * N)

        a = fwdpy11.Multiplicative(2.0)
        p = {
            "nregions": [],
            "sregions": [fwdpy11.ExpS(0, 1, 1, 0.01)],
            "recregions": [fwdpy11.Region(0, 1, 1)],
            "rates": (0.0, 0.00005, r),
            "gvalue": a,
            "prune_selected": True,
            "demography": fwdpy11.DiscreteDemography(),
            "simlen": 2 * N,
        }
        params = fwdpy11.ModelParams(**p)
        rng = fwdpy11.GSLrng(101 * 45 * 110 * 210)
        pop = fwdpy11.DiploidPopulation(N, 1.0)
        fwdpy11.evolvets(rng, pop, params, 100, track_mutation_counts=True)
        mc = fwdpy11.count_mutations(
            pop.tables, pop.mutations, [i for i in range(2 * pop.N)]
        )
        assert len(pop.fixations) > 0, "Test is meaningless without fixations"
        fixations = np.where(mc == 2 * pop.N)[0]
        self.assertEqual(len(fixations), 0)

        # Brute-force calculation of fixations
        brute_force = np.zeros(len(pop.mutations), dtype=np.int32)
        for g in pop.haploid_genomes:
            if g.n > 0:
                for k in g.smutations:
                    brute_force[k] += g.n

        self.assertTrue(np.array_equal(brute_force, mc))
        self.assertTrue(np.array_equal(brute_force, pop.mcounts))


class TestMetaData(unittest.TestCase):
    """
    Integration test.

    Rebuild genetic values from mutations
    stored in tree sequence
    """

    def testQtraitSim(self):
        N = 1000
        rho = 1.0
        r = rho / (4 * N)

        GSS = fwdpy11.GSS(VS=1.0, optimum=1.0)
        a = fwdpy11.Additive(2.0, GSS)
        p = {
            "nregions": [],
            "sregions": [fwdpy11.GaussianS(0, 1, 1, 0.25)],
            "recregions": [fwdpy11.PoissonInterval(0, 1, r)],
            "rates": (0.0, 0.005, None),
            "gvalue": a,
            "prune_selected": False,
            "demography": fwdpy11.DiscreteDemography(),
            "simlen": N,
        }
        params = fwdpy11.ModelParams(**p)
        rng = fwdpy11.GSLrng(101 * 45 * 110 * 210)
        pop = fwdpy11.DiploidPopulation(N, 1.0)

        # NOTE: we are trying to track the fate of C++
        # types based on their Python instances. We can
        # run into trouble as the underlying data are
        # over-written.  Thus, we store some of it
        # in a new class for testing.
        MD = namedtuple("MD", ["g", "e", "w", "label"])

        class Recorder(object):
            """ Records entire pop every 100 generations """

            def __init__(self):
                self.data = []

            def __call__(self, pop, recorder):
                if pop.generation % 100 == 0.0:
                    for i in pop.diploid_metadata:
                        self.data.append((pop.generation, MD(i.g, i.e, i.w, i.label)))
                    recorder.assign(np.arange(pop.N, dtype=np.int32))

        recorder = Recorder()
        fwdpy11.evolvets(rng, pop, params, 100, recorder)

        ancient_sample_metadata = np.array(pop.ancient_sample_metadata, copy=False)
        alive_sample_metadata = np.array(pop.diploid_metadata, copy=False)
        metadata = np.hstack((ancient_sample_metadata, alive_sample_metadata))

        nodes = np.array(pop.tables.nodes, copy=False)
        metadata_nodes = metadata["nodes"].flatten()
        metadata_node_times = nodes["time"][metadata_nodes]
        metadata_record_times = nodes["time"][metadata["nodes"][:, 0]]

        genetic_trait_values_from_sim = []
        genetic_values_from_ts = []
        for u in np.unique(metadata_node_times):
            samples_at_time_u = metadata_nodes[np.where(metadata_node_times == u)]
            vi = fwdpy11.VariantIterator(pop.tables, samples_at_time_u)
            sum_esizes = np.zeros(len(samples_at_time_u))
            for variant in vi:
                g = variant.genotypes
                r = variant.records[0]
                mutant = np.where(g == 1)[0]
                sum_esizes[mutant] += pop.mutations[r.key].s
            ind = int(len(samples_at_time_u) / 2)
            temp_gvalues = np.zeros(ind)
            temp_gvalues += sum_esizes[0::2]
            temp_gvalues += sum_esizes[1::2]
            genetic_values_from_ts.extend(temp_gvalues.tolist())
            genetic_trait_values_from_sim.extend(
                metadata["g"][np.where(metadata_record_times == u)[0]].tolist()
            )

        for i, j in zip(genetic_trait_values_from_sim, genetic_values_from_ts):
            self.assertAlmostEqual(i, j)

        # This is a second test that we correctly stored metadata
        # in the correct order
        for i, j in zip(pop.ancient_sample_metadata, recorder.data):
            self.assertEqual(i.g, j[1].g)
            self.assertEqual(i.e, j[1].e)
            self.assertEqual(i.w, j[1].w)
            self.assertEqual(i.label, j[1].label)


class TestDataMatrixIterator(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        # TODO add neutral variants
        self.N = 1000
        self.rho = 1.0
        self.theta = 100.0
        self.nreps = 500
        self.mu = self.theta / (4 * self.N)
        self.r = self.rho / (4 * self.N)

        self.GSS = fwdpy11.GSS(VS=1.0, optimum=0.0)
        a = fwdpy11.Additive(2.0, self.GSS)
        self.p = {
            "nregions": [],
            "sregions": [fwdpy11.GaussianS(0, 1, 1, 0.25)],
            "recregions": [fwdpy11.Region(0, 1, 1)],
            "rates": (0.0, 0.025, self.r),
            "gvalue": a,
            "prune_selected": False,
            "demography": fwdpy11.DiscreteDemography(),
            "simlen": self.N,
        }
        self.params = fwdpy11.ModelParams(**self.p)
        self.rng = fwdpy11.GSLrng(101 * 45 * 110 * 210)
        self.pop = fwdpy11.DiploidPopulation(self.N, 1.0)
        self.all_samples = [i for i in range(2 * self.N)]
        fwdpy11.evolvets(self.rng, self.pop, self.params, 100)
        self.dm = fwdpy11.data_matrix_from_tables(
            self.pop.tables, self.all_samples, True, True
        )
        self.neutral = np.array(self.dm.neutral)
        self.npos = np.array(self.dm.neutral.positions)
        self.selected = np.array(self.dm.selected)
        self.spos = np.array(self.dm.selected.positions)

    def test_entire_matrix(self):
        dmi = fwdpy11.DataMatrixIterator(
            self.pop.tables, self.all_samples, [(0, 1)], True, True
        )
        niterations = 0
        for dm in dmi:
            niterations += 1
            for i in dm.selected_keys:
                self.assertFalse(self.pop.mutations[i].neutral)
            self.assertTrue(np.array_equal(self.dm.selected_keys, dm.selected_keys))
            self.assertTrue(np.array_equal(dm.neutral, self.neutral))
            self.assertTrue(np.array_equal(dm.selected, self.selected))
        self.assertEqual(niterations, 1)

    def test_single_slice(self):
        dm = fwdpy11.DataMatrixIterator(
            self.pop.tables, self.all_samples, [(0.1, 0.2)], True, True
        )
        dm = next(dm)
        rows = np.where((self.spos >= 0.1) & (self.spos < 0.2))[0]
        pos_slice = self.spos[rows]
        selected_slice = self.selected[
            rows,
        ]
        self.assertTrue(np.array_equal(dm.selected_positions, pos_slice))
        self.assertTrue(np.array_equal(dm.selected, selected_slice))

    def test_nonoverlapping_slices(self):
        slices = [(0.1, 0.2), (0.21, 0.37), (0.5, 0.55)]
        dmi = fwdpy11.DataMatrixIterator(
            self.pop.tables, self.all_samples, slices, True, True
        )
        niterations = 0
        for r, dm in zip(slices, dmi):
            niterations += 1
            rows = np.where((self.spos >= r[0]) & (self.spos < r[1]))[0]
            pos_slice = self.spos[rows]
            selected_slice = self.selected[
                rows,
            ]
            self.assertTrue(np.array_equal(dm.selected_positions, pos_slice))
            self.assertTrue(np.array_equal(dm.selected, selected_slice))
        self.assertEqual(niterations, len(slices))

    def test_complex_slices(self):
        slices = [(0.1, 0.2), (0.15, 0.23), (0.21, 0.37), (0.38, 0.5337), (0.5, 0.55)]
        dmi = fwdpy11.DataMatrixIterator(
            self.pop.tables, self.all_samples, slices, True, True
        )
        niterations = 0
        for r, dm in zip(slices, dmi):
            niterations += 1
            rows = np.where((self.spos >= r[0]) & (self.spos < r[1]))[0]
            pos_slice = self.spos[rows]
            selected_slice = self.selected[
                rows,
            ]
            self.assertTrue(np.array_equal(dm.selected_positions, pos_slice))
            self.assertTrue(np.array_equal(dm.selected, selected_slice))
        self.assertEqual(niterations, len(slices))

    def test_nested_slices(self):
        slices = [
            (0.1, 0.2),
            (0.15, 0.19),
            (0.21, 0.37),
            (0.38, 0.5337),
            (0.39, 0.432),
            (0.5, 0.55),
        ]
        dmi = fwdpy11.DataMatrixIterator(
            self.pop.tables, self.all_samples, slices, True, True
        )
        for r, dm in zip(slices, dmi):
            rows = np.where((self.spos >= r[0]) & (self.spos < r[1]))[0]
            pos_slice = self.spos[rows]
            selected_slice = self.selected[
                rows,
            ]
            self.assertTrue(np.array_equal(dm.selected_positions, pos_slice))
            self.assertTrue(np.array_equal(dm.selected, selected_slice))


class TestTreeSequenceResettingDuringTimeSeriesAnalysis(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        class CountSamplesPerTimePoint(object):
            def __init__(self):
                self.sample_timepoints = []
                self.sample_sizes = []
                self.timepoint_seen = {}

            def __call__(self, pop):
                assert len(pop.tables.preserved_nodes) // 2 == len(
                    pop.ancient_sample_metadata
                )
                # Get the most recent ancient samples
                # and record their number.  We do this
                # by a "brute-force" approach
                for t, n, m in pop.sample_timepoints(False):
                    if t not in self.timepoint_seen:
                        self.timepoint_seen[t] = 1
                    else:
                        self.timepoint_seen[t] += 1
                    if t not in self.sample_timepoints:
                        self.sample_timepoints.append(t)
                        self.sample_sizes.append(len(n) // 2)

                    # simplify to each time point
                    tables, idmap = fwdpy11.simplify_tables(pop.tables, n)
                    for ni in n:
                        assert idmap[ni] != fwdpy11.NULL_NODE
                        assert tables.nodes[idmap[ni]].time == t

        self.N = 1000
        self.rho = 1.0
        self.r = self.rho / (4 * self.N)

        self.GSS = fwdpy11.GSS(VS=1.0, optimum=0.0)
        a = fwdpy11.Additive(2.0, self.GSS)
        self.p = {
            "nregions": [],
            "sregions": [fwdpy11.GaussianS(0, 1, 1, 0.25)],
            "recregions": [fwdpy11.Region(0, 1, 1)],
            "rates": (0.0, 0.025, self.r),
            "gvalue": a,
            "prune_selected": False,
            "demography": fwdpy11.DiscreteDemography(),
            "simlen": 101,
        }
        self.params = fwdpy11.ModelParams(**self.p)
        self.rng = fwdpy11.GSLrng(101 * 45 * 110 * 210)
        self.pop = fwdpy11.DiploidPopulation(self.N, 1.0)
        self.all_samples = [i for i in range(2 * self.N)]
        self.ancient_sample_recorder = fwdpy11.RandomAncientSamples(
            seed=42, samplesize=10, timepoints=[i for i in range(1, 101)]
        )
        self.resetter = CountSamplesPerTimePoint()
        fwdpy11.evolvets(
            self.rng,
            self.pop,
            self.params,
            5,
            recorder=self.ancient_sample_recorder,
            post_simplification_recorder=self.resetter,
        )

    def test_no_preserved_nodes(self):
        self.assertEqual(len(self.pop.tables.preserved_nodes), 0)

    def test_no_ancient_sample_metadata(self):
        self.assertEqual(len(self.pop.ancient_sample_metadata), 0)

    def test_all_timepoints_seen_exactly_once(self):
        for i, j in self.resetter.timepoint_seen.items():
            self.assertEqual(j, 1)

    def test_all_expected_timepoints_are_present(self):
        self.assertEqual(self.resetter.sample_timepoints, [i for i in range(1, 101)])

    def test_all_sample_sizes_are_correct(self):
        self.assertTrue(all([i == 10 for i in self.resetter.sample_sizes]) is True)


class TestDiploidPopulationInitialization(unittest.TestCase):
    def test_initial_nodes(self):
        pop = fwdpy11.DiploidPopulation(100, 1.0)
        md = np.array(pop.diploid_metadata, copy=False)
        n = md["nodes"].flatten()
        self.assertTrue(np.array_equal(n, np.arange(2 * pop.N, dtype=n.dtype)))


class TestTrackMutationCounts(unittest.TestCase):
    """
    Test for GitHub issue 326.

    Simplifying every generation w/strong positive selection is
    key to triggering the bug.
    """

    @classmethod
    def setUp(self):
        class FreqTracker(object):
            def __init__(self):
                self.freqs = dict()

            def __call__(self, pop, recorder):
                for i, j in enumerate(pop.mcounts):
                    if j == 2 * pop.N:
                        k = pop.mutations[i].key
                        if k in self.freqs:
                            if self.freqs[k][-1][1] < 2 * pop.N:
                                self.freqs[k].append((pop.generation, j))
                        else:
                            self.freqs[k] = [(pop.generation, j)]

        self.params, self.rng, self.pop = set_up_standard_pop_gen_model()
        self.ft = FreqTracker()

    def test_without_table_indexing(self):
        fwdpy11.evolvets(
            self.rng,
            self.pop,
            self.params,
            1,
            self.ft,
            track_mutation_counts=True,
            suppress_table_indexing=True,
        )
        assert len(self.pop.fixations) > 0, "Nothing fixed, so test case is not helpful"

        for i, j in zip(self.pop.fixations, self.pop.fixation_times):
            self.assertTrue(i.key in self.ft.freqs)
            self.assertEqual(self.ft.freqs[i.key][-1][0], j)

    def test_with_table_indexing(self):
        fwdpy11.evolvets(
            self.rng, self.pop, self.params, 1, self.ft, track_mutation_counts=True
        )
        assert len(self.pop.fixations) > 0, "Nothing fixed, so test case is not helpful"

        for i, j in zip(self.pop.fixations, self.pop.fixation_times):
            self.assertTrue(i.key in self.ft.freqs)
            self.assertEqual(self.ft.freqs[i.key][-1][0], j)


class TestRecapitation(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.params, self.rng, self.pop = set_up_standard_pop_gen_model(0.50)
        fwdpy11.evolvets(
            self.rng, self.pop, self.params, 100, preserve_first_generation=True
        )
        self.amd = np.array(self.pop.ancient_sample_metadata, copy=False)
        self.tskit_ts = self.pop.dump_tables_to_tskit()

    def test_ancient_sample_records(self):
        self.assertEqual(len(self.pop.ancient_sample_metadata), self.pop.N)
        self.assertEqual(len(self.pop.tables.preserved_nodes), 2 * self.pop.N)

    def test_ancient_sample_times(self):
        nodes = np.array(self.pop.tables.nodes, copy=False)
        self.assertTrue(np.all(nodes["time"][self.amd["nodes"].flatten()] == 0))

    def test_number_of_roots(self):
        tv = fwdpy11.TreeIterator(
            self.pop.tables,
            np.concatenate((self.pop.alive_nodes, self.amd["nodes"].flatten())),
        )
        for t in tv:
            self.assertEqual(len(t.roots), 2 * len(self.amd))

        for t in self.tskit_ts.trees():
            self.assertEqual(t.num_roots, 2 * len(self.amd))

    def test_dump_to_tskit(self):
        dumped_ts = self.pop.dump_tables_to_tskit()
        for i, j in enumerate(self.pop.ancient_sample_metadata):
            self.assertTrue(
                (
                    dumped_ts.tables.individuals.flags[
                        i + len(self.pop.diploid_metadata)
                    ]
                    & fwdpy11.tskit_tools.INDIVIDUAL_IS_PRESERVED
                )
                > 0
            )
            self.assertTrue(
                (
                    dumped_ts.tables.individuals.flags[
                        i + len(self.pop.diploid_metadata)
                    ]
                    & fwdpy11.tskit_tools.INDIVIDUAL_IS_FIRST_GENERATION
                )
                > 0
            )
            self.assertTrue(
                (
                    dumped_ts.tables.individuals.flags[
                        i + len(self.pop.diploid_metadata)
                    ]
                    & fwdpy11.tskit_tools.INDIVIDUAL_IS_ALIVE
                )
                == 0
            )
            for n in j.nodes:
                self.assertTrue(
                    (dumped_ts.tables.nodes.flags[n] & tskit.NODE_IS_SAMPLE) > 0
                )

    def test_recapitation(self):
        """
        Recapitation w/o recombination means we can predict
        exactly how many new nodes there will be.
        """
        coalesced_ts = msprime.simulate(
            Ne=self.pop.N,
            from_ts=self.tskit_ts,
            recombination_rate=self.params.rates.recombination_rate,
        )
        self.assertTrue(
            len(coalesced_ts.tables.nodes) > len(self.tskit_ts.tables.nodes)
        )

    def test_invalid_preservation(self):
        """
        Pop is already evolved, so we'll get an error now.
        """
        with self.assertRaises(ValueError):
            fwdpy11.evolvets(
                self.rng, self.pop, self.params, 100, preserve_first_generation=True
            )


class TestInvalidAttemptAtRecapitation(unittest.TestCase):
    def test_after_starting_from_msprime(self):
        ts = msprime.simulate(100, Ne=100)
        pop = fwdpy11.DiploidPopulation.create_from_tskit(ts)

        a = fwdpy11.Multiplicative(2.0)
        p = {
            "nregions": [],
            "sregions": [],
            "recregions": [],
            "rates": (0.0, 0.0, 0),
            "gvalue": a,
            "prune_selected": True,
            "demography": fwdpy11.DiscreteDemography(),
            "simlen": 100,
        }
        params = fwdpy11.ModelParams(**p)
        rng = fwdpy11.GSLrng(666 ** 2)
        with self.assertRaises(ValueError):
            fwdpy11.evolvets(rng, pop, params, 100, preserve_first_generation=True)


if __name__ == "__main__":
    unittest.main()
