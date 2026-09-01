#! /bin/csh -f

set test_dir = automated_testing/UPSY/component_test_PETSc_DMPLEX

rm -rf ${test_dir}/results*

mpiexec  -n 2  build/src/UPSY/UPSY_component_test_program_PETSc_DMPLEX  automated_testing/test_meshes_and_grids  ${test_dir}/results
find ${test_dir}/results -type f -name '*.h5' -exec sh -c 'python3 "$1" "$2" --output "${2%.h5}.png"' sh ${test_dir}/visualise_DMPLEX.py {} \;


