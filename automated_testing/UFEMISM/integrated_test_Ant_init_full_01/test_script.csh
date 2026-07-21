#! /bin/csh -f

set test_dir = automated_testing/UFEMISM/integrated_test_Ant_init_full_01

rm -rf ${test_dir}/results*

mpiexec  -n 2  UFEMISM_program  ${test_dir}/config.cfg
