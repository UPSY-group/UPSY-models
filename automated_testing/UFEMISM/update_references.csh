#!/bin/csh

# Define colors
set RED = "\033[0;31m"
set GREEN = "\033[0;32m"
set BLUE = "\033[0;34m"
set RESET = "\033[0m"

set all_tests = ( \
  integrated_test_SSA_icestream_small \
  integrated_test_ISMIP_HOM_small \
  integrated_test_Halfar_dome_small \
  integrated_test_MISMIP_mod_small \
  integrated_test_MISMIPplus_small \
  integrated_test_Berends2023nudging_exp1_small \
  integrated_test_Berends2023nudging_exp2_small \
  integrated_test_Ant_init_small_01 \
  integrated_test_ISMIP7_demo_small)

echo "Starting UFEMISM tests"

foreach test_name ($all_tests)

  printf "${GREEN}Starting UFEMISM test ${BLUE}${test_name}${GREEN}... ${RESET}\n"

  csh automated_testing/UFEMISM/${test_name}/test_script.csh
  rm -rf automated_testing/UFEMISM/${test_name}/reference
  mv automated_testing/UFEMISM/${test_name}/results_checksum   automated_testing/UFEMISM/${test_name}/reference
  find automated_testing/UFEMISM/${test_name}/results -type f -name '*.txt' -exec cp {} automated_testing/UFEMISM/${test_name}/reference \;

  if ($status != 0) then
      printf "${RED}Error: UFEMISM test ${BLUE}${test_name}${RED} failed with status $status ${RESET}\n"
      exit 1
  else
      printf "${GREEN}UFEMISM test ${BLUE}${test_name}${GREEN} completed successfully.${RESET}\n"
  endif
end

