#! /bin/csh -f

rm -rf automated_testing/UFEMISM/integrated_test_MISMIP_mod_small/results
mpiexec  -n 2 UFEMISM_program  automated_testing/UFEMISM/integrated_test_MISMIP_mod_small/config.cfg