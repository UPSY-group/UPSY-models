function(sanitize_fortran_sources output_var target_name)

  # CMake/Ninja must see generated files as real build outputs so it can track
  # timestamps on the original sources and regenerate sanitized copies only
  # when a source file actually changes.
  set(sanitized_sources "")
  get_property(existing_sanitized_sources GLOBAL PROPERTY SANITIZED_FORTRAN_OUTPUTS)
  if(NOT existing_sanitized_sources)
    set(existing_sanitized_sources "")
  endif()

  foreach(source_file IN LISTS ARGN)
    get_filename_component(source_file_abs "${source_file}" ABSOLUTE BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
    file(RELATIVE_PATH source_file_rel "${PROJECT_SOURCE_DIR}" "${source_file_abs}")
    set(sanitized_source "${PROJECT_BINARY_DIR}/sanitized_fortran/${target_name}/${source_file_rel}")

    list(FIND existing_sanitized_sources "${sanitized_source}" sanitized_source_index)
    if(sanitized_source_index EQUAL -1)
      add_custom_command(
        OUTPUT "${sanitized_source}"
        COMMAND "${CMAKE_COMMAND}"
        "-DINPUT_FILE:FILEPATH=${source_file_abs}"
        "-DOUTPUT_FILE:FILEPATH=${sanitized_source}"
        -P "${PROJECT_SOURCE_DIR}/cmake/sanitize_fortran_comment_quotes.cmake"
        DEPENDS "${source_file_abs}" "${PROJECT_SOURCE_DIR}/cmake/sanitize_fortran_comment_quotes.cmake"
        # The sanitizer removes quote characters from comments before the Fortran
        # compiler reads the source, which avoids cpp-style quote warnings.
        COMMENT "Sanitizing Fortran comments in ${source_file_rel}"
        VERBATIM
      )

      list(APPEND existing_sanitized_sources "${sanitized_source}")
      set_property(GLOBAL PROPERTY SANITIZED_FORTRAN_OUTPUTS "${existing_sanitized_sources}")
    endif()

    set_source_files_properties("${sanitized_source}" PROPERTIES GENERATED TRUE)

    list(APPEND sanitized_sources "${sanitized_source}")
  endforeach()

  set(${output_var} "${sanitized_sources}" PARENT_SCOPE)

endfunction()