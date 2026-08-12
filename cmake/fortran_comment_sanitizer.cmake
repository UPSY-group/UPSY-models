function(sanitize_fortran_sources output_var target_name)

  # CMake/Ninja must see generated files as real build outputs so it can track
  # timestamps on the original sources and regenerate sanitized copies only
  # when a source file actually changes.
  set(sanitized_sources "")

  foreach(source_file IN LISTS ARGN)
    get_filename_component(source_file_abs "${source_file}" ABSOLUTE BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
    file(RELATIVE_PATH source_file_rel "${PROJECT_SOURCE_DIR}" "${source_file_abs}")
    set(sanitized_source "${PROJECT_BINARY_DIR}/sanitized_fortran/${target_name}/${source_file_rel}")
    get_filename_component(sanitized_source_dir "${sanitized_source}" DIRECTORY)

    add_custom_command(
      OUTPUT "${sanitized_source}"
      COMMAND "${CMAKE_COMMAND}" -DINPUT_FILE=${source_file_abs} -DOUTPUT_FILE=${sanitized_source} -P "${PROJECT_SOURCE_DIR}/cmake/sanitize_fortran_comment_apostrophes.cmake"
      DEPENDS "${source_file_abs}" "${PROJECT_SOURCE_DIR}/cmake/sanitize_fortran_comment_apostrophes.cmake"
      # The sanitizer removes quote characters from comments before the Fortran
      # compiler reads the source, which avoids cpp-style quote warnings.
      COMMENT "Sanitizing Fortran comments in ${source_file_rel}"
      VERBATIM
    )

    list(APPEND sanitized_sources "${sanitized_source}")
  endforeach()

  set(${output_var} "${sanitized_sources}" PARENT_SCOPE)

endfunction()