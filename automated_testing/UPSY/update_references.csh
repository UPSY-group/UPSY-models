#!/bin/csh

set all_tests = ( \
  component_test_mesh_creation \
  component_test_mesh_discretisation \
  component_test_mesh_focussing \
  component_test_mesh_remapping_mesh_grid \
  component_test_mesh_remapping_mesh_mesh)

foreach test_name ($all_tests)

  csh automated_testing/UPSY/${test_name}/test_script.csh
  rm -rf automated_testing/UPSY/${test_name}/reference
  mv automated_testing/UPSY/${test_name}/results_checksum   automated_testing/UPSY/${test_name}/reference
  cp automated_testing/UPSY/${test_name}/results/*.txt      automated_testing/UPSY/${test_name}/reference

end