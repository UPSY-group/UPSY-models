#! /bin/csh -f

set test_dir = automated_testing/UPSY/component_test_remapping

# Start from scratch: the weights are always recalculated (rather than read from
# UPSY_remapping_weights/), so that the test covers UPSY_calculate_remapping_weights
rm -rf ${test_dir}/results ${test_dir}/results_checksum ${test_dir}/weights ${test_dir}/figures

# Figures are only for visual inspection when running the test locally, so skip them on CI
if ($?CI) then
  python3 ${test_dir}/remapping_test.py  automated_testing/test_meshes_and_grids  ${test_dir}/results  ${test_dir}/weights
else
  python3 ${test_dir}/remapping_test.py  automated_testing/test_meshes_and_grids  ${test_dir}/results  ${test_dir}/weights  --figures ${test_dir}/figures
endif
if ($status != 0) exit 1

python3 automated_testing/reduce_all_netcdfs_in_folder_to_checksum.py ${test_dir}
