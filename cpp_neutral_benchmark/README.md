# Benchmarking program

This directory contains a C++-only program for simple benchmarking:

* No selection
* No mutation
* Recombination is allowed
* Tree sequence recording happens

This program exists because it is easier to use `perf` on it than on the entire Python package.

## Compiling

The build system is "GNU pedantic", so you need to specify everything:

```sh
automake --add-missing
autoreconf
CXXFLAGS="-O3 -g -DNDEBUG" ./configure
make
```

The boost program options library is needed at compile/run time.
