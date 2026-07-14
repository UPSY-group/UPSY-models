#! /bin/csh -f

set test_dir = automated_testing/UFEMISM/integrated_test_Halfar_dome_full

rm -rf ${test_dir}/results*

mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_Halfar_40km.cfg
mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_Halfar_20km.cfg
mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_Halfar_10km.cfg
mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_Halfar_5km.cfg

mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_Halfar_static_40km.cfg
mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_Halfar_static_20km.cfg
mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_Halfar_static_10km.cfg
mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_Halfar_static_5km.cfg
