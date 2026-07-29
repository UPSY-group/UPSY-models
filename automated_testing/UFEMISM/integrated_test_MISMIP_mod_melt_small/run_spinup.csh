#! /bin/csh -f

set test_dir = automated_testing/UFEMISM/integrated_test_MISMIP_mod_melt_small

rm -rf ${test_dir}/spun_up_geometry
mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_spinup.cfg
