#!/usr/bin/env bash

set -euo pipefail

mpiexec  -n 2  \
  build/src/UPSY/UPSY_unit_test_program             automated_testing/UPSY/unit_tests/results
mpiexec  -n 7  --map-by :OVERSUBSCRIBE  \
  build/src/UPSY/UPSY_multinode_unit_test_program   automated_testing/UPSY/unit_tests_multinode/results
csh ./automated_testing/check_unit_tests.csh
