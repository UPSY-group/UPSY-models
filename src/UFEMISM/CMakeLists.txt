# Collect LADDIE/UFEMISM source files
file( GLOB_RECURSE LADDIE_sources_raw  "${PROJECT_SOURCE_DIR}/src/LADDIE/*.f90")
list(REMOVE_ITEM LADDIE_sources_raw "${PROJECT_SOURCE_DIR}/src/LADDIE/main/LADDIE_program.f90")

file( GLOB_RECURSE UFEMISM_sources_raw "${CMAKE_CURRENT_SOURCE_DIR}/*.f90")
file( GLOB UFEMISM_program_sources "${CMAKE_CURRENT_SOURCE_DIR}/main/*_program.f90")

foreach(program_source IN LISTS UFEMISM_program_sources)
	list(REMOVE_ITEM UFEMISM_sources_raw "${program_source}")
endforeach()

# Create static library for UFEMISM/LADDIE sources
add_library( UFEMISM_static_library STATIC ${LADDIE_sources_raw} ${UFEMISM_sources_raw})
target_link_libraries( UFEMISM_static_library PUBLIC UPSY_static_library)
link_to_dependencies( UFEMISM_static_library)
add_compiler_flags( UFEMISM_static_library)

target_include_directories( UFEMISM_static_library PUBLIC ${PROJECT_SOURCE_DIR}/include ${CMAKE_CURRENT_BINARY_DIR}/mod_files/UFEMISM_static_library)
set_property( TARGET UFEMISM_static_library PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_files/UFEMISM_static_library)

# Initialise an empty list of executables
set( all_programs "")

# UFEMISM programs (main and tests)
foreach(program_source IN LISTS UFEMISM_program_sources)
	get_filename_component( program_name "${program_source}" NAME_WE)
	add_executable( ${program_name} "${program_source}")
	list( APPEND all_programs ${program_name})
endforeach()

# Link all executables to the UFEMISM static library and external dependencies,
# and set compiler flags and dedicated Fortran module output folders.
foreach(prog IN LISTS all_programs)
	target_link_libraries(${prog} PUBLIC UFEMISM_static_library)
	link_to_dependencies(${prog})
	add_compiler_flags(${prog})
	set_property(TARGET ${prog} PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_files/${prog})
endforeach()

# Print a green summary line once all UFEMISM programs built successfully.
if(all_programs)
	string(ASCII 27 ESC)
	add_custom_target( UFEMISM_all_programs_build_success ALL
		DEPENDS ${all_programs}
		COMMAND ${CMAKE_COMMAND} -E echo "${ESC}[32mUFEMISM build success: all UFEMISM programs compiled and linked.${ESC}[0m"
		COMMENT "Checking UFEMISM program build status"
	)
endif()
