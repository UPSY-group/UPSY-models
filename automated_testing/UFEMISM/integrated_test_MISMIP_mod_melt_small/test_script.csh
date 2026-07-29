#! /bin/csh -f

set test_dir = automated_testing/UFEMISM/integrated_test_MISMIP_mod_melt_small

# rm -rf ${test_dir}/results*

# mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_spinup_part0_40km.cfg
# mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_spinup_10km.cfg
# mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_melt_40km.cfg
# mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_melt_10km.cfg

rm -rf ${test_dir}/results_melt_10km

mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_melt_10km.cfg
