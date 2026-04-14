#! /bin/csh -f

mpiexec  -n 2  build/src/UPSY/UPSY_component_test_program_mesh_creation

rm -rf automated_testing/UPSY/component_test_mesh_creation/results
mv  test_meshes_and_grids  automated_testing/UPSY/component_test_mesh_creation/results