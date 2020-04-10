#
# Copyright (C) 2019 Kevin Thornton <krthornt@uci.edu>
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
import datetime
import sys
from collections import namedtuple

import numpy as np
import pandas as pd
import seaborn as sns

import fwdpy11

THETA = 0.5e3
NSAM_PER_DEME = 10
N_PER_DEME = 1000

Datum = namedtuple("Datum", ["migrate", "Fst"])


def make_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--nreps",
        type=int,
        default=100,
        help="Number of replicates to simulate per migration rate",
    )
    parser.add_argument(
        "--filename", type=str, default="migtest.png", help="File name for final plot"
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=666,
        help="Random number seed, used for both fwdpy11 and for numpy",
    )
    parser.add_argument("--quickrun", default=False, action="store_true")
    return parser


def plot_results(df, filename, nreps):
    """
    Jitter plot of the raw data + mean and expected values
    """
    # df.sort_values(by=['migrate'])

    df["migrate_str"] = ["{:.2E}".format(i) for i in df["migrate"]]
    df.sort_values(by=["migrate"])

    g = df.groupby(["migrate"]).mean().reset_index()
    g["efst"] = 1.0 / (1 + 4 * 2 * N_PER_DEME * g["migrate"])
    g["migrate_str"] = ["{:.2E}".format(i) for i in g["migrate"]]

    order = [(i, j) for i, j in zip(g["migrate"], g["migrate_str"])]
    order = sorted(order, key=lambda x: x[0])
    order = [i[1] for i in order]

    fourNm = [(i, "{:.2E}".format(4 * N_PER_DEME * i)) for i in g["migrate"]]
    fourNm = sorted(fourNm, key=lambda x: x[0])
    fourNm = [i[1] for i in fourNm]

    ax = sns.stripplot(
        x="migrate_str", y="Fst", alpha=0.5, data=df, jitter=True, order=order
    )
    now = datetime.datetime.now()
    day = now.strftime("%D")
    title = "Run on {}, based on {} replicates per migration rate.".format(day, nreps)
    ax.set_ylim(0, 1)
    ax.set_xlabel(r"$4Nm$")
    ax.set_ylabel(r"$F_{st}$")
    ax.plot("migrate_str", "Fst", data=g, label=r"mean $F_{st}$")
    ax.plot("migrate_str", "efst", data=g, label=r"expected $F_{st}$")
    ax.set_xticklabels(fourNm)
    ax.legend()
    ax.set_title(title, fontdict={"fontsize": 10})
    ax.get_figure().savefig(filename)


def fst(m):
    """
    Lazy implementation of Hudson, Slatkin, and
    Maddison's Fst.

    Hudson, R. R., M. Slatkin, and W. P. Maddison. 1992.
    “Estimation of Levels of Gene Flow from DNA Sequence Data.”
    Genetics 132 (2): 583–89.
    """
    nc = 2 * NSAM_PER_DEME

    # pi w/in each deme
    dcounts0 = np.sum(m[:, :nc], axis=1)
    dcounts1 = np.sum(m[:, nc:], axis=1)
    pi0 = np.sum(dcounts0 * (nc - dcounts0))
    pi0 *= 2.0 / (nc * (nc - 1))
    pi1 = np.sum(dcounts1 * (nc - dcounts1))
    pi1 *= 2.0 / (nc * (nc - 1))

    # pi in total sample
    dcounts = np.sum(m, axis=1)
    piT = np.sum(dcounts * (m.shape[1] - dcounts))
    piT *= 2.0 / (m.shape[1] * (m.shape[1] - 1))

    # Fst
    piS = (pi0 + pi1) / 2.0
    piD = 2.0 * (piT - piS)
    return piD / (piD + piS)


def runsim(seed, mij):
    migmatrix = np.identity(2)
    migmatrix[:] = mij
    np.fill_diagonal(migmatrix, 0)
    rs = np.sum(migmatrix, axis=1)
    np.fill_diagonal(migmatrix, 1.0 - rs)
    np.random.seed(seed)
    rng = fwdpy11.GSLrng(seed)
    pop = fwdpy11.DiploidPopulation(2 * N_PER_DEME, 1.0)
    mm = [fwdpy11.move_individuals(0, 0, 1, 0.5)]

    dd = fwdpy11.DiscreteDemography(
        mass_migrations=mm, migmatrix=fwdpy11.MigrationMatrix(migmatrix)
    )

    pdict = {
        "sregions": [],
        "nregions": [],
        "recregions": [fwdpy11.PoissonInterval(0, 1, 1e2 / (4 * N_PER_DEME))],
        "gvalue": fwdpy11.Multiplicative(2.0),
        "rates": (0, 0, None),
        "demography": dd,
        "simlen": 20 * N_PER_DEME,
    }

    params = fwdpy11.ModelParams(**pdict)

    fwdpy11.evolvets(rng, pop, params, 100, suppress_table_indexing=True)
    fwdpy11.infinite_sites(rng, pop, THETA / float(4 * N_PER_DEME))
    md = np.array(pop.diploid_metadata, copy=False)
    demes = np.unique(md["deme"])
    w = np.where(md["deme"] == demes[0])[0]
    s0 = np.random.choice(w, NSAM_PER_DEME, replace=False)
    w = np.where(md["deme"] == demes[1])[0]
    s1 = np.random.choice(w, NSAM_PER_DEME, replace=False)
    samples = np.concatenate((md["nodes"][s0].flatten(), md["nodes"][s1].flatten()))
    dm = fwdpy11.data_matrix_from_tables(pop.tables, samples, True, False)
    ngm = np.array(dm.neutral, copy=False)
    return Datum(mij, fst(ngm))


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])
    if args.quickrun is True:
        N_PER_DEME = 200
        max_workers = 1
    else:
        max_workers = None

    rng = fwdpy11.GSLrng(args.seed)
    np.random.seed(args.seed)
    migrates = np.linspace(0, 2.0 / (4 * N_PER_DEME), 6)
    migrates[0] = 1e-5
    mij = [*migrates.tolist()] * args.nreps
    seeds = np.random.randint(
        0, np.iinfo(np.uint32).max, len(mij) * args.nreps, dtype=np.uint32
    )
    results = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as pool:
        futures = {pool.submit(runsim, i, j) for i, j in zip(seeds, mij)}
        for fut in concurrent.futures.as_completed(futures):
            result = fut.result()
            results.append(result)

    df = pd.DataFrame(results, columns=Datum._fields)
    plot_results(df, args.filename, args.nreps)
