#! /bin/csh -f

set test_dir = automated_testing/UFEMISM/integrated_test_ISMIP_HOM_full
set config_dir = ${test_dir}/all_ISMIP_HOM_config_files

rm -rf ${test_dir}/results*
mkdir  ${test_dir}/results

set experiments    = (A B C D)
set length_scales  = (160 80 40 20 10 5)
set Stokes_approxs = (SIASSA DIVA BPA)

rm -rf $test_dir/results_ISMIP_HOM*

foreach experiment ($experiments)
  foreach length_scale ($length_scales)
    foreach Stokes_approx ($Stokes_approxs)

      set exp_name = ${experiment}_${length_scale}_${Stokes_approx}
      set exp_output_dir = $test_dir/results_ISMIP_HOM_${exp_name}

      rm -rf $exp_output_dir

      mpiexec -n 2 UFEMISM_program $config_dir/config_ISMIP_HOM_${exp_name}.cfg

      mv $exp_output_dir/transect_ISMIP-HOM.nc ${test_dir}/results/transect_${exp_name}.nc

      rm -rf $exp_output_dir

    end
  end
end
