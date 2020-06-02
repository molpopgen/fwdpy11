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
"""
Local adaptation of a quantitative trait to differing optima.
"""

import argparse
import math
import sys
from collections import namedtuple

import numpy as np
import pandas as pd

import fwdpy11

# Simulations with tree sequence recording need
# to know the max position in a genome.  Here,
# we use a length of 1.0. Thus, all mutation
# and recombination events will be uniform
# random variables on the continuous interval
# [0, GENOME_LENGTH).
GENOME_LENGTH = 1.0

# When recording quant-genetic statistics during a simulation,
# we will use this type. Named tuples are extremely efficient,
# and they are easily converted into Pandas DataFrame objects,
# which is very convenient for analysis and output.
Datum = namedtuple("Data", ["generation", "deme", "gbar", "vg", "wbar"])


def make_parser():
    """
    Create a command-line interface to the script.
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    required = parser.add_argument_group("Required arguments")
    required.add_argument("--popsize", "-N", type=int, help="Diploid population size")
    required.add_argument(
        "--mu", "-m", type=float, help="Mutation rate (per gamete, per generation)"
    )
    required.add_argument(
        "--sigma",
        "-s",
        type=float,
        help="Standard deviation of Gaussian" "distribution of mutational effects",
    )

    optional = parser.add_argument_group("Optional arguments")
    optional.add_argument(
        "--rho", type=float, default=1000.0, help="Scaled recombination rate, rho=4Nr"
    )
    optional.add_argument(
        "--VS",
        type=float,
        default=10.0,
        help="Inverse strength of stabilizing selection",
    )
    optional.add_argument(
        "--opt", type=float, default=1.0, help="Value of new phenotypic optimum"
    )
    optional.add_argument(
        "--migrates",
        type=float,
        nargs=2,
        default=None,
        help="Migration rates from 0 to 1 and 1 to 0, respectively.",
    )
    optional.add_argument(
        "--time",
        type=float,
        default=0.1,
        help="Amount of time to simulate past" "optimum shift, in units of N",
    )
    optional.add_argument(
        "--plotfile", type=str, default=None, help="File name for plot"
    )
    optional.add_argument("--seed", type=int, default=42, help="Random number seed.")

    return parser


def validate_arguments(args):
    """
    Validate input arguments.
    Note: this is likely incomplete.
    """
    if args.popsize is None:
        raise ValueError("popsize cannot be None")
    if args.mu < 0:
        raise ValueError("mu must be non-negative")
    if args.mu is None:
        raise ValueError("mu cannot be None")
    if args.mu < 0 or math.isfinite(args.mu) is False:
        raise ValueError("Mutation rate must be non-negative and finite")
    if args.sigma is None:
        raise ValueError("sigma cannot be none")
    if args.sigma < 0 or math.isfinite(args.sigma) is False:
        raise ValueError(
            "Std. dev. of distribution of effect sizes"
            "must be non-negative and finite"
        )

    if args.migrates is not None:
        for m in args.migrates:
            if m < 0 or m > 1:
                raise ValueError("migration rates must be 0 <= m <= 1")


def make_migmatrix(migrates):
    if migrates is None:
        return None
    mm = np.zeros(4).reshape(2, 2)
    mm[0, 1] = migrates[1]
    mm[1, 0] = migrates[0]
    rs = np.sum(mm, axis=1)
    np.fill_diagonal(mm, 1.0 - rs)
    return fwdpy11.MigrationMatrix(mm)


class Recorder(object):
    """
    fwdpy11 allows you to define objects that record data
    from populations during simulation.  Such objects must
    be callable, and the easiest way to do things is to
    create a class with a __call__ function.
    """

    def __init__(self, start):
        self.data = []
        self.start = start

    def __call__(self, pop, recorder):
        if pop.generation >= self.start:
            # Record mean trait value each generation.
            md = np.array(pop.diploid_metadata, copy=False)
            demes = np.unique(md["deme"])
            for d in demes:
                w = np.where(md["deme"] == d)[0]
                gbar = md["g"][w].mean()
                vg = md["g"][w].var()
                wbar = md["w"][w].mean()
                self.data.append(Datum(pop.generation, d, gbar, vg, wbar))


def plot_output(data, filename):
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    fig = plt.figure(figsize=(9, 3))
    gs = gridspec.GridSpec(ncols=3, nrows=1, figure=fig)
    ax_gbar = fig.add_subplot(gs[0, 0])
    ax_vg = fig.add_subplot(gs[0, 1])
    ax_wbar = fig.add_subplot(gs[0, 2])

    df = pd.DataFrame(data, columns=Datum._fields)
    g = df.groupby(["deme"])
    for n, gi in g:
        ax_gbar.plot(gi["generation"], gi["gbar"], label="Deme {}".format(n))
        ax_vg.plot(gi["generation"], gi["vg"], label="Deme {}".format(n))
        ax_wbar.plot(gi["generation"], gi["wbar"], label="Deme {}".format(n))

    for ax in [ax_gbar, ax_vg, ax_wbar]:
        ax.set_xlabel("Generation")

    ax_gbar.set_ylabel(r"$\bar{g}$")
    ax_vg.set_ylabel(r"$V(G)$")
    ax_wbar.set_ylabel(r"$\bar{w}$")
    ax_gbar.legend()
    plt.tight_layout()
    plt.savefig(filename)


def runsim(args):
    """
    Run the simulation.
    """
    pop = fwdpy11.DiploidPopulation(2 * args.popsize, GENOME_LENGTH)

    np.random.seed(args.seed)
    rng = fwdpy11.GSLrng(args.seed)

    GSSmo0 = fwdpy11.GSSmo(
        [
            fwdpy11.Optimum(when=0, optimum=0.0, VS=args.VS),
            fwdpy11.Optimum(when=10 * args.popsize, optimum=args.opt, VS=args.VS),
        ]
    )
    GSSmo1 = fwdpy11.GSSmo(
        [
            fwdpy11.Optimum(when=0, optimum=0.0, VS=args.VS),
            fwdpy11.Optimum(
                when=10 * args.popsize, optimum=-1.0 * args.opt, VS=args.VS
            ),
        ]
    )

    mm = make_migmatrix(args.migrates)
    dd = fwdpy11.DiscreteDemography(
        mass_migrations=[fwdpy11.move_individuals(0, 0, 1, 0.5)], migmatrix=mm
    )

    p = {
        "nregions": [],  # No neutral mutations -- add them later!
        "gvalue": [fwdpy11.Additive(2.0, GSSmo0), fwdpy11.Additive(2.0, GSSmo1)],
        "sregions": [fwdpy11.GaussianS(0, GENOME_LENGTH, 1, args.sigma)],
        "recregions": [fwdpy11.Region(0, GENOME_LENGTH, 1)],
        "rates": (0.0, args.mu, args.rho / float(4 * args.popsize)),
        # Keep mutations at frequency 1 in the pop if they affect fitness.
        "prune_selected": False,
        "demography": dd,
        "simlen": 10 * args.popsize + int(args.popsize * args.time),
    }
    params = fwdpy11.ModelParams(**p)

    r = Recorder(10 * args.popsize)
    fwdpy11.evolvets(rng, pop, params, 100, r, suppress_table_indexing=True)
    if args.plotfile is not None:
        plot_output(r.data, args.plotfile)


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])
    validate_arguments(args)
    runsim(args)
