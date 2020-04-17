import os

import ycm_core

cwd = os.path.dirname(os.path.realpath(__file__))

FLAGS = [
    "-Wall",
    "-Weffc++",
    "-std=c++14",
    "-Werr",
    "-x",
    "c++",
    "-I",
    ".",
    "-I",
    cwd + "/fwdpy11/headers/fwdpp",
    "-I",
    cwd + "/fwdpy11/headers",
]


def FlagsForFile(filename):
    dirname = os.path.dirname(filename)
    flags = {"flags": FLAGS, "do_cache": True}
    return flags
