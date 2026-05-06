#! /bin/csh -f

set test_dir = automated_testing/UFEMISM/integrated_test_Berends2023nudging_exp1_small

rm -rf ${test_dir}/results*

mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_01_exp_I_spinup_40km_part0.cfg
mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_02_exp_I_spinup_40km.cfg
mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_03_exp_I_spinup_10km.cfg
mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_05_exp_I_inversion_10km_H_dHdt_flowline.cfg
mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_06_exp_I_inversion_10km_H_dHdt_local.cfg
mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_07_exp_I_inversion_10km_H_u_flowline.cfg

mkdir ${test_dir}/results

set subtests = ( \
  results_01_exp_I_spinup_40km_part0 \
  results_02_exp_I_spinup_40km \
  results_03_exp_I_spinup_10km \
  results_05_exp_I_inversion_5km_H_dHdt_flowline \
  results_06_exp_I_inversion_5km_H_dHdt_local \
  results_07_exp_I_inversion_5km_H_u_flowline)

foreach subtest ($subtests)
  mv ${test_dir}/${subtest}/main_output_ANT_00001.nc   ${test_dir}/results/main_output_ANT_${subtest}_00001.nc
  mv ${test_dir}/${subtest}/main_output_ANT_grid.nc    ${test_dir}/results/main_output_ANT_${subtest}_grid.nc
  mv ${test_dir}/${subtest}/checksum_logfile.txt       ${test_dir}/results/checksum_logfile_${subtest}.txt

  rm -rf ${test_dir}/${subtest}
end

.venv/bin/python3  automated_testing/reduce_all_netcdfs_in_folder_to_checksum.py  ${test_dir}