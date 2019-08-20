"""
Evolution of a single genomic region under Gaussian Stabilizing Selection
with a moving optimum, or "GSSmo".

The simulation allows the recording of quantitative-genetics statistics
over time as well as the preservation of samples into the tree sequence
after the optimum shift.

The demographic model is constant-size Wright-Fisher for 10N generations
around an optimal trait value of zero.  Then, the optimum changes and
the population continues to evolve for a length of time specified
by the user.

The mutation model is uniform across the genome with effect sizes
given by a Gaussian distribution with mean zero and a user-specified
standard deviation.

The genetic model is additive effects with fitness determined by
Gaussian stabilizing selection based on the squared distance from
the optimum.
"""

import fwdpy11
from collections import namedtuple
import pandas as pd
import sqlite3
import math
import numpy as np
import pickle
import lzma
import sys
import os
import argparse

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
DATA = namedtuple('Data', ['generation', 'zbar', 'vg', 'wbar'])


def make_parser():
    """
    Create a command-line interface to the script.
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group("Required arguments")
    required.add_argument("--popsize", "-N", type=int,
                          help="Diploid population size")
    required.add_argument("--mu", "-m", type=float,
                          help="Mutation rate (per gamete, per generation)")
    required.add_argument("--sigma", "-s", type=float,
                          help="Standard deviation of Gaussian"
                          "distribution of mutational effects")
    required.add_argument("--filename", "-f", type=str, default=None,
                          help="Output file name.  The population will be"
                          "pickled to this file, and the file will be"
                          "compressed using the lzma algorithm")

    optional = parser.add_argument_group("Optional arguments")
    optional.add_argument("--rho", type=float, default=1000.0,
                          help="Scaled recombination rate, rho=4Nr")
    optional.add_argument("--VS", type=float, default=1.0,
                          help="Inverse strength of stabilizing selection")
    optional.add_argument("--opt", type=float, default=1.0,
                          help="Value of new phenotypic optimum")
    optional.add_argument("--time", type=float, default=0.1,
                          help="Amount of time to simulate past"
                          "optimum shift, in units of N")
    optional.add_argument("--record", "-r", action="store_true",
                          help="Record statistics such as VG from population "
                          "each generation")
    optional.add_argument('--statfile', '-S', type=str, default=None,
                          help="File name to record statistics."
                          "The format is an sqlite3 database.")
    optional.add_argument("--preserve", type=int, default=-1,
                          help="Record ancient samples every X generations"
                          "after optimum shift. A value of -1 means not"
                          "to record.")
    optional.add_argument("--num_ind", '-n', type=int, default=50,
                          help="Number of diploids to record as"
                          "ancient samples")
    optional.add_argument("--seed", type=int, default=42,
                          help="Random number seed.")

    return parser


def validate_arguments(args):
    """
    Validate input arguments.
    Note: this is likely incomplete.
    """
    if args.popsize < 0:
        raise ValueError("Population size must be non-negative")
    if args.mu < 0 or math.isfinite(args.mu) is False:
        raise ValueError("Mutation rate must be non-negative and finite")
    if args.sigma < 0 or math.isfinite(args.sigma) is False:
        raise ValueError("Std. dev. of distribution of effect sizes"
                         "must be non-negative and finite")
    if args.filename is None:
        raise ValueError("ouptut filename cannot be None")

    if args.record > 0 and args.statfile is None:
        raise ValueError("Must profile a stats file name when recording"
                         "statistics during simulation")

    if args.statfile is not None and os.path.exists(args.statfile):
        raise RuntimeError("statfile already exists!")

    if args.preserve > 0:
        if args.num_ind > args.popsize:
            raise ValueError("Number of ancient samples is"
                             "greater than population size")


class Recorder(object):
    """
    fwdpy11 allows you to define objects that record data
    from populations during simulation.  Such objects must
    be callable, and the easiest way to do things is to
    create a class with a __call__ function.
    """

    def __init__(self, record_stats, interval, nsam):
        self.data = []
        self.record = record_stats
        self.interval = interval
        self.nsam = nsam

    def __call__(self, pop, recorder):
        if self.record is True:
            # Record mean trait value each generation.
            t = np.array(pop.diploid_metadata, copy=False)
            self.data.append(DATA(pop.generation, t['g'].mean(),
                                  t['g'].var(), t['w'].mean()))

        if self.interval > 0 and pop.generation >= 10*pop.N:
            if pop.generation % self.interval == 0.0:
                if self.nsam < pop.N:
                    s = np.random.choice(pop.N, self.nsam, replace=False)
                else:
                    s = np.arange(pop.N, dtype=np.int32)

                recorder.assign(s)


def runsim(args):
    """
    Run the simulation and deliver output to files.
    """
    pop = fwdpy11.DiploidPopulation(args.popsize, GENOME_LENGTH)

    rng = fwdpy11.GSLrng(args.seed)

    GSSmo = fwdpy11.GSSmo(
        [(0, 0, args.VS), (10*args.popsize, args.opt, args.VS)])

    popsizes = [args.popsize]*(10*args.popsize +
                               int(args.time*float(args.popsize)))
    p = {'nregions': [],  # No neutral mutations -- add them later!
         'gvalue': fwdpy11.Additive(2.0, GSSmo),
         'sregions': [fwdpy11.GaussianS(0, GENOME_LENGTH, 1, args.sigma)],
         'recregions': [fwdpy11.Region(0, GENOME_LENGTH, 1)],
         'rates': (0.0, args.mu, args.rho/float(4*args.popsize)),
         # Keep mutations at frequency 1 in the pop if they affect fitness.
         'prune_selected': False,
         'demography':  np.array(popsizes, dtype=np.uint32)
         }
    params = fwdpy11.ModelParams(**p)

    r = Recorder(args.record, args.preserve, args.num_ind)
    fwdpy11.evolvets(rng, pop, params, 100, r, suppress_table_indexing=True)

    with lzma.open(args.filename, 'wb') as f:
        pickle.dump(pop, f)

    if args.statfile is not None:
        stats = pd.DataFrame(r.data, columns=DATA._fields)
        # Write the statistics to an sqlite3 database,
        # which can be processed in R via dplyr.
        conn = sqlite3.connect(args.statfile)
        stats.to_sql('data', conn)


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])
    validate_arguments(args)
    runsim(args)
