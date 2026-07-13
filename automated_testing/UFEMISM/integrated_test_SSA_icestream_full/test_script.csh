#! /bin/csh -f

set test_dir = automated_testing/UFEMISM/integrated_test_SSA_icestream_full

rm -rf ${test_dir}/results*

mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_01_32km.cfg
mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_02_16km.cfg
mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_03_8km.cfg
mpiexec  -n 2  UFEMISM_program  ${test_dir}/config_04_4km.cfg
