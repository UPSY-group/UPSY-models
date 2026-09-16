function(link_to_dependencies target)

  # Detect platform
  if(APPLE)
      set(IS_MACOS TRUE)
  elseif(UNIX)
      set(IS_LINUX TRUE)
  endif()

  # =============
  # == OpenMPI ==
  # =============

  # Find MPI package
  find_package(MPI REQUIRED Fortran)

  # Add include directories and link libraries
  target_link_libraries(${target} PRIVATE MPI::MPI_Fortran)

  # ===========
  # == PETSc ==
  # ===========

  # Try pkg-config first: this covers both conda environments (which set
  # PKG_CONFIG_PATH on activation) and HPC module systems that also export a
  # PETSc.pc (e.g. EasyBuild-based stacks such as Snellius).
  find_package(PkgConfig)
  set(PETSC_FOUND_VIA_PKGCONFIG FALSE)
  if(PKG_CONFIG_FOUND)
    pkg_check_modules(PETSC QUIET PETSc)
    if(PETSC_FOUND)
      set(PETSC_FOUND_VIA_PKGCONFIG TRUE)
    endif()
  endif()

  if(PETSC_FOUND_VIA_PKGCONFIG)

    include_directories(${PETSC_INCLUDE_DIRS})
    link_directories(${PETSC_LIBRARY_DIRS})
    add_definitions(${PETSC_CFLAGS_OTHER})

    if(IS_LINUX)
      find_package(HDF5 REQUIRED COMPONENTS C HL Fortran)
      target_link_libraries(${target} PRIVATE ${HDF5_LIBRARIES})
      include_directories(${HDF5_INCLUDE_DIRS})
      target_link_libraries(${target} PRIVATE ${PETSC_LIBRARIES})
    elseif(IS_MACOS)
        target_link_libraries(${target} PRIVATE ${PETSC_LIBRARY_DIRS}/libpetsc.dylib)
    endif()

  else()

    # No PETSc.pc available (e.g. no conda environment, and the loaded PETSc
    # module doesn't export one). Fall back to the PETSC_DIR (and optional
    # PETSC_ARCH) environment variables that `module load PETSc` sets on HPC
    # systems such as Snellius.
    if(NOT DEFINED ENV{PETSC_DIR})
      message(FATAL_ERROR "Could not find PETSc: no PETSc.pc found via pkg-config, and the "
                           "PETSC_DIR environment variable is not set. Either activate the conda "
                           "environment, or load the PETSc module (e.g. `module load PETSc`).")
    endif()

    set(PETSC_DIR_ENV "$ENV{PETSC_DIR}")
    set(PETSC_HPC_INCLUDE_DIRS "${PETSC_DIR_ENV}/include")
    set(PETSC_HPC_LIBRARY_DIRS "${PETSC_DIR_ENV}/lib")
    if(DEFINED ENV{PETSC_ARCH} AND NOT "$ENV{PETSC_ARCH}" STREQUAL "")
      list(APPEND PETSC_HPC_INCLUDE_DIRS "${PETSC_DIR_ENV}/$ENV{PETSC_ARCH}/include")
      list(APPEND PETSC_HPC_LIBRARY_DIRS "${PETSC_DIR_ENV}/$ENV{PETSC_ARCH}/lib")
    endif()

    include_directories(${PETSC_HPC_INCLUDE_DIRS})
    link_directories(${PETSC_HPC_LIBRARY_DIRS})
    target_link_libraries(${target} PRIVATE petsc)

  endif()

  # ============
  # == NetCDF ==
  # ============

  find_package(PkgConfig)
  set(NETCDF_FOUND_VIA_PKGCONFIG FALSE)
  if(PKG_CONFIG_FOUND)
    pkg_check_modules(NETCDF QUIET netcdf-fortran)
    if(NETCDF_FOUND)
      set(NETCDF_FOUND_VIA_PKGCONFIG TRUE)
    endif()
  endif()

  if(NETCDF_FOUND_VIA_PKGCONFIG)

    include_directories(${NETCDF_INCLUDE_DIRS})
    link_directories(${NETCDF_LIBRARY_DIRS})
    add_definitions(${NETCDF_CFLAGS_OTHER})

    if(IS_LINUX)
        target_link_libraries(${target} PRIVATE ${NETCDF_LIBRARIES})
    elseif(IS_MACOS)
        target_link_libraries(${target} PRIVATE ${NETCDF_LIBRARY_DIRS}/libnetcdff.dylib)
    endif()

  else()

    # No netcdf-fortran.pc available. Fall back to nf-config, which ships
    # with every NetCDF-Fortran install (including HPC modules, e.g. the
    # `netCDF-Fortran` module on Snellius) and reports the correct
    # include/link flags directly.
    find_program(NF_CONFIG_EXECUTABLE nf-config)
    if(NOT NF_CONFIG_EXECUTABLE)
      message(FATAL_ERROR "Could not find NetCDF-Fortran: no netcdf-fortran.pc found via "
                           "pkg-config, and nf-config is not on PATH. Either activate the conda "
                           "environment, or load the NetCDF-Fortran module "
                           "(e.g. `module load netCDF-Fortran`).")
    endif()

    execute_process(COMMAND ${NF_CONFIG_EXECUTABLE} --includedir
                     OUTPUT_VARIABLE NETCDF_HPC_INCLUDE_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${NF_CONFIG_EXECUTABLE} --flibs
                     OUTPUT_VARIABLE NETCDF_HPC_FLIBS OUTPUT_STRIP_TRAILING_WHITESPACE)

    include_directories(${NETCDF_HPC_INCLUDE_DIR})
    separate_arguments(NETCDF_HPC_FLIBS_LIST UNIX_COMMAND "${NETCDF_HPC_FLIBS}")
    target_link_libraries(${target} PRIVATE ${NETCDF_HPC_FLIBS_LIST})

  endif()

endfunction()