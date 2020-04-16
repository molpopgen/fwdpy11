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
import argparse
import concurrent.futures
import sys
from collections import namedtuple

import msprime
import numpy as np
import pandas as pd

import fwdpy11


SimOutput = namedtuple("SimOutput", ["S2N", "Pi2N", "Sn", "Pin"])


def make_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("-N", type=int, help="Diploid population size")
    parser.add_argument("--rho", type=float, help="4Nr")
    parser.add_argument("--theta", type=float, help="4Nu")
    parser.add_argument(
        "--simlen", type=int, help="Generations to run the forward simulation"
    )
    parser.add_argument("--seed", type=int, help="Random number seed")
    parser.add_argument("--nreps", type=int, help="Number of simulation replicates")
    parser.add_argument("--nsam", type=int, help="Sample size to analyze")

    optional = parser.add_argument_group("Optional arguments")
    optional.add_argument("--model", type=str, help="msprime model", default="dtwf")
    return parser


def run_msprime(Ne, rho, theta, model, seed, input_ts=None):
    if input_ts is None:
        ts = msprime.simulate(
            2 * Ne,
            Ne=Ne,
            model=model,
            random_seed=seed,
            recombination_rate=rho / Ne / 4,
            mutation_rate=theta / Ne / 4,
        )
    else:
        ts = msprime.simulate(
            from_ts=input_ts,
            Ne=Ne,
            model=model,
            random_seed=seed,
            recombination_rate=rho / Ne / 4,
        )

    return ts


def run_msprime_get_stats(Ne, nsam, rho, theta, model, seed):
    ts = run_msprime(Ne, rho, theta, model, seed)
    samples = [i for i in ts.samples()]
    S1 = len(ts.tables.sites)
    pi1 = ts.diversity([samples])
    rsample = np.random.choice(samples, nsam, replace=False)
    S2 = ts.segregating_sites([rsample])
    pi2 = ts.diversity([rsample])
    return SimOutput(S1, pi1, S2, pi2)


def run_sim(args):
    user_args, msprime_seed, fwdpy11_seed = args
    np.random.seed(msprime_seed)
    pop = fwdpy11.DiploidPopulation(user_args.N, 1.0)

    pdict = {
        "nregions": [],
        "sregions": [],
        "recregions": [
            fwdpy11.PoissonInterval(
                0, pop.tables.genome_length, user_args.rho / pop.N / 4
            )
        ],
        "rates": (0, 0, None),
        "gvalue": fwdpy11.Multiplicative(1),
        "demography": fwdpy11.DiscreteDemography(),
        "simlen": user_args.simlen,
    }
    params = fwdpy11.ModelParams(**pdict)
    rng = fwdpy11.GSLrng(fwdpy11_seed)
    fwdpy11.evolvets(rng, pop, params, 100, preserve_first_generation=True)
    ts = pop.dump_tables_to_tskit()
    completed_ts = run_msprime(
        user_args.N, user_args.rho, 0.0, user_args.model, msprime_seed, ts
    )
    samples = pop.alive_nodes
    completed_ts = msprime.mutate(
        completed_ts, rate=user_args.theta / user_args.N / 4, random_seed=msprime_seed
    )
    S2N = completed_ts.segregating_sites([samples])[0]
    pi2N = completed_ts.diversity([samples])[0]
    rsample = np.random.choice(samples, user_args.nsam, replace=False)
    Sn = completed_ts.segregating_sites([rsample])[0]
    pin = completed_ts.diversity([rsample])[0]

    return SimOutput(S2N, pi2N, Sn, pin)


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])

    np.random.seed(args.seed)

    msprime_seeds = []
    fwdpy11_seeds = []

    for i in range(args.nreps):
        candidate = np.random.randint(0, np.iinfo(np.uint32).max)
        while candidate in msprime_seeds:
            candidate = np.random.randint(0, np.iinfo(np.uint32).max)
        msprime_seeds.append(candidate)
        candidate = np.random.randint(0, np.iinfo(np.uint32).max)
        while candidate in fwdpy11_seeds:
            candidate = np.random.randint(0, np.iinfo(np.uint32).max)
        fwdpy11_seeds.append(candidate)

    results = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = {
            executor.submit(run_sim, (args, i, j))
            for i, j in zip(msprime_seeds, fwdpy11_seeds)
        }
        for future in concurrent.futures.as_completed(futures):
            results.append(future.result())

    results = pd.DataFrame(results, columns=SimOutput._fields)
    results["watterson_n"] = results.Sn / (1.0 / (np.arange(args.nsam - 1) + 1)).sum()
    print(f"Means from fwdpy11:\n{results.mean()}")

    msprime_seeds = []

    for i in range(args.nreps):
        candidate = np.random.randint(0, np.iinfo(np.uint32).max)
        while candidate in msprime_seeds:
            candidate = np.random.randint(0, np.iinfo(np.uint32).max)
        msprime_seeds.append(candidate)
    msprime_results = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = {
            executor.submit(
                run_msprime_get_stats,
                args.N,
                args.nsam,
                args.rho,
                args.theta,
                args.model,
                i,
            )
            for i in msprime_seeds
        }
        for future in concurrent.futures.as_completed(futures):
            msprime_results.append(future.result())
    msprime_results = pd.DataFrame(msprime_results, columns=SimOutput._fields)
    msprime_results["watterson_n"] = (
        msprime_results.Sn / (1.0 / (np.arange(args.nsam - 1) + 1)).sum()
    )
    print(f"Means from msprime:\n{msprime_results.mean()}")
