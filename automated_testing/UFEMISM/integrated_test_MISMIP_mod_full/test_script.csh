#! /bin/csh -f

set test_dir = automated_testing/UFEMISM/integrated_test_MISMIP_mod_full

rm -rf ${test_dir}/results*

mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_01_spinup_40km.cfg
mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_02_spinup_10km.cfg
mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_03_advance_10km.cfg
mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_04_retreat_10km.cfg
