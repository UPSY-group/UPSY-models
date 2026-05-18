#!/bin/csh

# Define colors
set RED = "\033[0;31m"
set GREEN = "\033[0;32m"
set BLUE = "\033[0;34m"
set RESET = "\033[0m"

set all_tests = ( \
  component_test_mesh_creation \
  component_test_mesh_discretisation \
  component_test_mesh_focussing \
  component_test_mesh_remapping_mesh_grid \
  component_test_mesh_remapping_mesh_mesh)

foreach test_name ($all_tests)

  printf "${GREEN}Starting UPSY test ${BLUE}${test_name}${GREEN}... ${RESET}\n"

  csh automated_testing/UPSY/${test_name}/test_script.csh
  rm -rf automated_testing/UPSY/${test_name}/reference
  mv automated_testing/UPSY/${test_name}/results_checksum   automated_testing/UPSY/${test_name}/reference
  find automated_testing/UPSY/${test_name}/results -type f -name '*.txt' -exec cp {} automated_testing/UPSY/${test_name}/reference \;

  if ($status != 0) then
      printf "${RED}Error: UPSY test ${BLUE}${test_name}${RED} failed with status $status ${RESET}\n"
      exit 1
  else
      printf "${GREEN}UPSY test ${BLUE}${test_name}${GREEN} completed successfully.${RESET}\n"
  endif
end
