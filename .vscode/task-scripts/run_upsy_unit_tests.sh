#!/usr/bin/env bash

set -euo pipefail

mpiexec -n 2 ./build/src/UPSY/UPSY_unit_test_program
mpiexec -n 7 ./build/src/UPSY/UPSY_multinode_unit_test_program
csh ./automated_testing/check_unit_tests.csh
