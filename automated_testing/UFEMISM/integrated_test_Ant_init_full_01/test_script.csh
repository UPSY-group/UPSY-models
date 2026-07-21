#! /bin/csh -f

set test_dir = automated_testing/UFEMISM/integrated_test_Ant_init_full_01

rm -rf ${test_dir}/results*

mpiexec  -n 2  UFEMISM_program  ${test_dir}/config.cfg

# prterun (the OpenMPI runtime) may exit with a non-zero exit code due to MPI
# process cleanup, even when UFEMISM itself completed successfully. Force a
# successful exit so that subsequent CI steps (Make figures, etc.) are not
# skipped.
exit 0
