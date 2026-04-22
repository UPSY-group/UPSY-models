#! /bin/csh -f

rm -rf automated_testing/UFEMISM/integrated_test_SSA_icestream_small/results
mpiexec  -n 2  UFEMISM_program  automated_testing/UFEMISM/integrated_test_SSA_icestream_small/config.cfg