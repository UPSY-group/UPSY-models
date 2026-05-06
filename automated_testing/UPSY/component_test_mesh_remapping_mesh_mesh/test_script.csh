#! /bin/csh -f

set test_dir = automated_testing/UPSY/component_test_mesh_remapping_mesh_mesh

mpiexec  -n 2  build/src/UPSY/UPSY_component_test_program_mesh_remapping_mesh_mesh  automated_testing/test_meshes_and_grids  ${test_dir}/results

.venv/bin/python3 automated_testing/reduce_all_netcdfs_in_folder_to_checksum.py ${test_dir}