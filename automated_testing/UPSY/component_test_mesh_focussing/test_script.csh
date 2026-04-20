#! /bin/csh -f

mpiexec  -n 2  build/src/UPSY/UPSY_component_test_program_mesh_focussing  automated_testing/test_meshes_and_grids  automated_testing/UPSY/component_test_mesh_focussing/results
