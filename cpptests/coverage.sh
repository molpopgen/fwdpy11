#!/bin/bash

gcov *.cc
lcov -c -d . -o coverage_unfiltered.info
lcov -e coverage_unfiltered.info -o coverage_remove_boost.info '*fwdpy11*'
lcov -r coverage_remove_boost.info -o coverage.info '*cpptests*'
genhtml coverage.info -o fwdpy11_coverage
