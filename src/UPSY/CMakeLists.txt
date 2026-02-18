# Collect UPSY/LADDIE/UFEMISM source files
file( GLOB_RECURSE UPSY_sources_raw    "${CMAKE_CURRENT_SOURCE_DIR}/*.f90")

# Create static library for UPSY
add_library( UPSY_static_library STATIC ${UPSY_sources_raw})
link_to_dependencies( UPSY_static_library)
add_compiler_flags( UPSY_static_library)

target_include_directories(UPSY_static_library PUBLIC ${PROJECT_SOURCE_DIR}/include ${CMAKE_CURRENT_BINARY_DIR}/mod_files)

set_property(TARGET UPSY_static_library PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_files)

# Initialise an empty list of executables
set( all_programs "")

# UPSY programs
add_executable( UPSY_unit_test_program           "${CMAKE_CURRENT_SOURCE_DIR}/UPSY_unit_test_program.f90")
add_executable( UPSY_multinode_unit_test_program "${CMAKE_CURRENT_SOURCE_DIR}/UPSY_multinode_unit_test_program.f90")
add_executable( UPSY_component_test_program      "${CMAKE_CURRENT_SOURCE_DIR}/UPSY_component_test_program.f90")
add_executable( UPSY_component_test_program_mesh_creation      "${CMAKE_CURRENT_SOURCE_DIR}/UPSY_component_test_program_mesh_creation.f90")
list( APPEND all_programs UPSY_unit_test_program)
list( APPEND all_programs UPSY_multinode_unit_test_program)
list( APPEND all_programs UPSY_component_test_program)
list( APPEND all_programs UPSY_component_test_program_mesh_creation)

# Link all executables to the UPSY static library and the external dependencies,
# and set the compiler flags
foreach(prog IN LISTS all_programs)
    target_link_libraries(${prog} PUBLIC UPSY_static_library)
    link_to_dependencies(${prog})
    add_compiler_flags(${prog})
endforeach()
