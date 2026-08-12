# Collect LADDIE/UFEMISM source files
file( GLOB_RECURSE LADDIE_sources_raw  "${CMAKE_CURRENT_SOURCE_DIR}/*.f90")
list(REMOVE_ITEM LADDIE_sources_raw "${CMAKE_CURRENT_SOURCE_DIR}/main/LADDIE_program.f90")
sanitize_fortran_sources(LADDIE_sources_sanitized LADDIE_program ${LADDIE_sources_raw})
sanitize_fortran_sources(LADDIE_program_source_sanitized LADDIE_program "${CMAKE_CURRENT_SOURCE_DIR}/main/LADDIE_program.f90")

file( GLOB_RECURSE UFEMISM_sources_raw "${PROJECT_SOURCE_DIR}/src/UFEMISM/*.f90")
file( GLOB UFEMISM_program_sources "${PROJECT_SOURCE_DIR}/src/UFEMISM/main/*_program.f90")

foreach(program_source IN LISTS UFEMISM_program_sources)
        list(REMOVE_ITEM UFEMISM_sources_raw "${program_source}")
endforeach()
sanitize_fortran_sources(UFEMISM_sources_sanitized LADDIE_program ${UFEMISM_sources_raw})

# LADDIE program
add_executable( LADDIE_program ${LADDIE_sources_sanitized} ${UFEMISM_sources_sanitized} ${LADDIE_program_source_sanitized})

target_link_libraries( LADDIE_program PUBLIC UPSY_static_library)
link_to_dependencies(LADDIE_program)
add_compiler_flags(LADDIE_program)
