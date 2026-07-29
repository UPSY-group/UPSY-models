#! /bin/csh -f

set test_dir = automated_testing/UFEMISM/integrated_test_MISMIP_mod_melt_small

# rm -rf ${test_dir}/results_melt
# mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_melt.cfg

rm -rf ${test_dir}/results_melt_hires
mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_melt_hires.cfg
