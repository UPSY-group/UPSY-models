#!/bin/csh

# Define colors
set RED = "\033[0;31m"
set GREEN = "\033[0;32m"
set BLUE = "\033[0;34m"
set RESET = "\033[0m"

set os_name = `uname -s`
if ("$os_name" == "Darwin") then
  set reference_dir = "reference_macos"
else if ("$os_name" == "Linux") then
  set reference_dir = "reference_linux"
else
  printf "${RED}Error: unsupported OS \"${os_name}\".${RESET}\n"
  exit 1
endif

# Use a direct solver (MUMPS LU) instead of an iterative Krylov solver, so that small
# floating-point differences cannot change the number of solver iterations (and thereby
# the solution). PETSc applies these options in KSPSetFromOptions, overriding the
# KSP/PC types from the config. Only applies to test runs, not to regular simulations.
# NOTE: not set for the UPSY tests, since component_test_PETSc_matrix_solving
#       specifically tests the different Krylov solvers and preconditioners.
setenv PETSC_OPTIONS "-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -ksp_initial_guess_nonzero false"

# Force single-threaded BLAS/OpenMP, so that results are reproducible across machines
# with different core counts (threaded BLAS reductions change the floating-point
# summation order). Only applies to test runs, not to regular simulations.
setenv OMP_NUM_THREADS 1
setenv OPENBLAS_NUM_THREADS 1

set all_tests = ( \
  integrated_test_SSA_icestream_small \
  integrated_test_SSA_notime_MISMIP_mod_full \
  integrated_test_ISMIP_HOM_small \
  integrated_test_Halfar_dome_small \
  integrated_test_MISMIP_mod_small \
  integrated_test_MISMIPplus_small \
  integrated_test_Berends2023nudging_exp1_small \
  integrated_test_Berends2023nudging_exp2_small \
  integrated_test_Ant_init_small_01 \
  integrated_test_ISMIP7_demo_small)

echo "Starting UFEMISM tests (updating ${reference_dir})"

foreach test_name ($all_tests)

  printf "${GREEN}Starting UFEMISM test ${BLUE}${test_name}${GREEN}... ${RESET}\n"

  csh automated_testing/UFEMISM/${test_name}/test_script.csh
  rm -rf automated_testing/UFEMISM/${test_name}/${reference_dir}
  mv automated_testing/UFEMISM/${test_name}/results_checksum   automated_testing/UFEMISM/${test_name}/${reference_dir}
  find automated_testing/UFEMISM/${test_name}/results -type f -name '*.txt' -exec cp {} automated_testing/UFEMISM/${test_name}/${reference_dir} \;

  if ($status != 0) then
      printf "${RED}Error: UFEMISM test ${BLUE}${test_name}${RED} failed with status $status ${RESET}\n"
      exit 1
  else
      printf "${GREEN}UFEMISM test ${BLUE}${test_name}${GREEN} completed successfully.${RESET}\n"
  endif
end

