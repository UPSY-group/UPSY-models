#!/bin/csh

set all_tests = ( \
  integrated_test_SSA_icestream_small \
  integrated_test_ISMIP_HOM_small \
  integrated_test_Halfar_dome_small \
  integrated_test_MISMIP_mod_small \
  integrated_test_MISMIPplus_small \
  integrated_test_Berends2023nudging_exp1_small \
  integrated_test_Berends2023nudging_exp2_small \
  integrated_test_Ant_init_small_01)

foreach test_name ($all_tests)

  csh automated_testing/UFEMISM/${test_name}/test_script.csh
  rm -rf automated_testing/UFEMISM/${test_name}/reference
  mv automated_testing/UFEMISM/${test_name}/results_checksum   automated_testing/UFEMISM/${test_name}/reference
  cp automated_testing/UFEMISM/${test_name}/results/*.txt      automated_testing/UFEMISM/${test_name}/reference

end