#! /bin/csh -f

set test_dir = automated_testing/UFEMISM/integrated_test_Halfar_dome_small

rm -rf ${test_dir}/results*

mpiexec  -n 2 UFEMISM_program  ${test_dir}/config_Halfar_40km.cfg
mpiexec  -n 2 UFEMISM_program  ${test_dir}/config_Halfar_static_40km.cfg

mkdir ${test_dir}/results

mv ${test_dir}/results_Halfar_40km/main_output_ANT_00001.nc          ${test_dir}/results/main_output_ANT_Halfar_40km_00001.nc
mv ${test_dir}/results_Halfar_40km/main_output_ANT_grid.nc           ${test_dir}/results/main_output_ANT_Halfar_40km_grid.nc
mv ${test_dir}/results_Halfar_40km/checksum_logfile.txt              ${test_dir}/results/checksum_logfile_Halfar_40km_grid.txt

rm -rf ${test_dir}/results_Halfar_40km

mv ${test_dir}/results_Halfar_static_40km/main_output_ANT_00001.nc   ${test_dir}/results/main_output_ANT_Halfar_static_40km_00001.nc
mv ${test_dir}/results_Halfar_static_40km/main_output_ANT_grid.nc    ${test_dir}/results/main_output_ANT_Halfar_static_40km_grid.nc
mv ${test_dir}/results_Halfar_static_40km/checksum_logfile.txt       ${test_dir}/results/checksum_logfile_Halfar_static_40km_grid.txt

rm -rf ${test_dir}/results_Halfar_static_40km

.venv/bin/python3 automated_testing/reduce_all_netcdfs_in_folder_to_checksum.py ${test_dir}