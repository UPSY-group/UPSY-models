#! /bin/csh -f

set test_dir = automated_testing/UPSY/component_test_ugrid

mpiexec  -n 2  build/src/UPSY/UPSY_component_test_program_ugrid  automated_testing/test_meshes_and_grids  ${test_dir}/results

rm -rf ${test_dir}/results_checksum
python3 automated_testing/reduce_all_netcdfs_in_folder_to_checksum.py ${test_dir}

# Validate the ugrid files and save plots of them
rm -rf ${test_dir}/figures
foreach filename (${test_dir}/results/ugrid_*.nc)
  python3 tools/python/validate_ugrid.py ${filename} --save-dir ${test_dir}/figures
end
