#! /bin/csh -f

mpiexec  -n 7  --map-by :OVERSUBSCRIBE  build/src/UPSY/UPSY_multinode_unit_test_program  automated_testing/UPSY/unit_tests_multinode/results
