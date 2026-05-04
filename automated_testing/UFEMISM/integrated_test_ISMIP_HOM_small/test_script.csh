#! /bin/csh -f

set test_dir = automated_testing/UFEMISM/integrated_test_ISMIP_HOM_small
set config_dir = ${test_dir}/all_ISMIP_HOM_config_files

rm -rf ${test_dir}/results*
mkdir  ${test_dir}/results

set experiments    = (A)                 # Full: (A B C D)             # Only do experiment A in this small test
set length_scales  = (160 80 40)         # Full: (160 80 40 20 10 5)   # Only do the first three length scales
set Stokes_approxs = (SIASSA DIVA BPA)   # Full: (SIASSA DIVA BPA)     # Do all three Stokes approximations

rm -rf $test_dir/results_ISMIP_HOM*

foreach experiment ($experiments)
  foreach length_scale ($length_scales)
    foreach Stokes_approx ($Stokes_approxs)

      set exp_name = ${experiment}_${length_scale}_${Stokes_approx}
      set exp_output_dir = $test_dir/results_ISMIP_HOM_${exp_name}

      rm -rf $exp_output_dir

      mpiexec -n 2 UFEMISM_program $config_dir/config_ISMIP_HOM_${exp_name}.cfg

      mv $exp_output_dir/main_output_ANT_00001.nc ${test_dir}/results/results_ISMIP_HOM_${exp_name}_mesh.nc
      mv $exp_output_dir/main_output_ANT_grid.nc  ${test_dir}/results/results_ISMIP_HOM_${exp_name}_grid.nc
      mv $exp_output_dir/checksum_logfile.txt     ${test_dir}/results/checksum_logfile_${exp_name}.txt

      rm -rf $exp_output_dir

    end
  end
end

.venv/bin/python3 automated_testing/reduce_all_netcdfs_in_folder_to_checksum.py ${test_dir}