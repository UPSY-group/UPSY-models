#! /bin/csh -f

set test_dir = automated_testing/UFEMISM/integrated_test_MISMIP_mod_melt_small

rm -rf ${test_dir}/results_melt
mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_melt.cfg

# python3 automated_testing/reduce_all_netcdfs_in_folder_to_checksum.py ${test_dir}
