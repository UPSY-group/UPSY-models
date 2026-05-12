#!/bin/csh

set all_tests = ( \
  component_test_mesh_creation \
  component_test_mesh_discretisation \
  component_test_mesh_focussing \
  component_test_mesh_remapping_mesh_grid \
  component_test_mesh_remapping_mesh_mesh)

foreach test_name ($all_tests)

  csh automated_testing/UPSY/${test_name}/test_script.csh
  echo "Finished script" ${test_name}
  rm -rf automated_testing/UPSY/${test_name}/reference
  echo "Removed reference" ${test_name}
  mv automated_testing/UPSY/${test_name}/results_checksum   automated_testing/UPSY/${test_name}/reference
  echo "Moved results_checksum" ${test_name}
  find automated_testing/UPSY/${test_name}/results -type f -name '*.txt' -exec cp {} automated_testing/UPSY/${test_name}/reference \;
  echo "Copied reference" ${test_name}

end
