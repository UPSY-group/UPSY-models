#! /bin/csh -f

set test_dir = automated_testing/UPSY/component_test_PETSc_matrix_solving

mpiexec  -n 2  build/src/UPSY/UPSY_component_test_program_PETSc_matrix_solving  ${test_dir}/results

python3 automated_testing/reduce_all_netcdfs_in_folder_to_checksum.py ${test_dir}
