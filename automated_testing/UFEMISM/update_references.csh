#!/bin/csh

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

  csh automated_testing/UFEMISM/${test_name}/test_script.csh
  echo "Finished script" ${test_name}
  rm -rf automated_testing/UFEMISM/${test_name}/reference
  echo "Removed reference" ${test_name}
  mv automated_testing/UFEMISM/${test_name}/results_checksum   automated_testing/UFEMISM/${test_name}/reference
  echo "Moved results_checksum" ${test_name}
  find automated_testing/UFEMISM/${test_name}/results -type f -name '*.txt' -exec cp {} automated_testing/UFEMISM/${test_name}/reference \;
  echo "Copied reference" ${test_name}

end
