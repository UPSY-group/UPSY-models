#! /bin/csh -f

set test_dir = automated_testing/UFEMISM/integrated_test_MISMIPplus_full

rm -rf ${test_dir}/results*

mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_01_5km_spinup_part0.cfg
mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_02_5km_spinup.cfg
mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_03_5km_ice1r.cfg
mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_04_4km_spinup.cfg
mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_05_4km_ice1r.cfg

#python3 automated_testing/reduce_all_netcdfs_in_folder_to_checksum.py ${test_dir}
