# Collect UPSY/LADDIE/UFEMISM source files
file( GLOB_RECURSE UPSY_sources_raw    "${CMAKE_CURRENT_SOURCE_DIR}/*.f90")
sanitize_fortran_sources(UPSY_sources_sanitized UPSY_static_library ${UPSY_sources_raw})

# Create static library for UPSY
add_library( UPSY_static_library STATIC ${UPSY_sources_sanitized})
link_to_dependencies( UPSY_static_library)
add_compiler_flags( UPSY_static_library)

if(EXISTS ${PROJECT_SOURCE_DIR}/include)
	target_include_directories(UPSY_static_library PUBLIC ${PROJECT_SOURCE_DIR}/include)
endif()
target_include_directories(UPSY_static_library PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/mod_files)

set_property(TARGET UPSY_static_library PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_files)

# Initialise an empty list of executables
set( all_programs "")

set( UPSY_program_sources_raw
	"${CMAKE_CURRENT_SOURCE_DIR}/UPSY_unit_test_program.f90"
	"${CMAKE_CURRENT_SOURCE_DIR}/UPSY_multinode_unit_test_program.f90"
	"${CMAKE_CURRENT_SOURCE_DIR}/UPSY_component_test_program_PETSc_matrix_solving.f90"
	"${CMAKE_CURRENT_SOURCE_DIR}/UPSY_component_test_program_mesh_creation.f90"
	"${CMAKE_CURRENT_SOURCE_DIR}/UPSY_component_test_program_mesh_focussing.f90"
	"${CMAKE_CURRENT_SOURCE_DIR}/UPSY_component_test_program_mesh_discretisation.f90"
	"${CMAKE_CURRENT_SOURCE_DIR}/UPSY_component_test_program_mesh_remapping_mesh_grid.f90"
	"${CMAKE_CURRENT_SOURCE_DIR}/UPSY_component_test_program_mesh_remapping_mesh_mesh.f90")

sanitize_fortran_sources(UPSY_program_sources_sanitized UPSY_programs ${UPSY_program_sources_raw})

# UPSY programs
foreach(program_source IN LISTS UPSY_program_sources_sanitized)
		get_filename_component(program_name "${program_source}" NAME_WE)
		add_executable(${program_name} "${program_source}")
		list(APPEND all_programs ${program_name})
endforeach()

# Link all executables to the UPSY static library and the external dependencies,
# and set the compiler flags
foreach(prog IN LISTS all_programs)
    target_link_libraries(${prog} PUBLIC UPSY_static_library)
    link_to_dependencies(${prog})
    add_compiler_flags(${prog})
endforeach()

# Print a green summary line once all UPSY programs built successfully.
if(all_programs)
	string(ASCII 27 ESC)
	add_custom_target( UPSY_all_programs_build_success ALL
		DEPENDS ${all_programs}
		COMMAND ${CMAKE_COMMAND} -E echo "${ESC}[32mUPSY build success: all UPSY programs compiled and linked.${ESC}[0m"
		COMMENT "Checking UPSY program build status"
	)
endif()
