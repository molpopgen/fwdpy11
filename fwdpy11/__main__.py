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

import argparse
import sys


def print_includes():
    from ._dev import get_fwdpp_includes
    from ._dev import get_includes

    idirs = [get_includes(), get_fwdpp_includes()]

    print(" ".join("-I" + idir for idir in idirs))


def get_fwdpy11_includes():
    from ._dev import get_includes

    return get_includes()


def get_fwdpp_includes():
    from ._dev import get_fwdpp_includes

    return get_fwdpp_includes()


def print_mako():
    from ._dev import minimal_mako

    print(minimal_mako())


def main():
    parser = argparse.ArgumentParser(prog="python -m fwdpy11")
    parser.add_argument(
        "--includes",
        action="store_true",
        help="Print the CPPFLAGS required to use fwdpy11 and fwdpp headers.",
    )
    parser.add_argument(
        "--fwdpy11_headers", action="store_true", help="Get location of fwdpy11 headers"
    )
    parser.add_argument(
        "--fwdpp_headers", action="store_true", help="Get location of fwdpp headers"
    )
    parser.add_argument(
        "--mako",
        action="store_true",
        help="Print minimal mako header for use with cppimport.",
    )
    args = parser.parse_args()
    if not sys.argv[1:]:
        parser.print_help()
    if args.includes:
        print_includes()
        return
    if args.fwdpy11_headers:
        print(get_fwdpy11_includes())
        return
    if args.fwdpp_headers:
        print(get_fwdpp_includes())
        return
    if args.mako:
        print_mako()
        return


if __name__ == "__main__":
    main()
