if(NOT DEFINED INPUT_FILE)
  message(FATAL_ERROR "INPUT_FILE is required")
endif()

if(NOT DEFINED OUTPUT_FILE)
  message(FATAL_ERROR "OUTPUT_FILE is required")
endif()

file(READ "${INPUT_FILE}" source_text)
string(REPLACE "\r\n" "\n" source_text "${source_text}")
string(REPLACE "\r" "\n" source_text "${source_text}")

# Walk the file line by line instead of converting it to a CMake list. Comment
# text can contain semicolons, and CMake would otherwise split those into extra
# list elements before we have a chance to sanitize them.
set(output_text "")

set(remaining_text "${source_text}")

while(NOT remaining_text STREQUAL "")
  string(FIND "${remaining_text}" "\n" newline_index)
  if(newline_index EQUAL -1)
    set(source_line "${remaining_text}")
    set(remaining_text "")
  else()
    string(SUBSTRING "${remaining_text}" 0 ${newline_index} source_line)
    math(EXPR next_index "${newline_index} + 1")
    string(LENGTH "${remaining_text}" remaining_length)
    if(next_index LESS remaining_length)
      string(SUBSTRING "${remaining_text}" ${next_index} -1 remaining_text)
    else()
      set(remaining_text "")
    endif()
  endif()

  if(source_line MATCHES "^[ ]*!")
    # For full-line comments, strip both kinds of quotes and keep the line
    # otherwise unchanged so the build output stays readable.
    string(REPLACE "'" "" processed_line "${source_line}")
    string(REPLACE "\"" "" processed_line "${processed_line}")
    string(APPEND output_text "${processed_line}\n")
  else()
    set(processed_line "")
    set(in_single_quote FALSE)
    set(in_comment FALSE)

    string(LENGTH "${source_line}" line_length)
    set(index 0)

    while(index LESS line_length)
      string(SUBSTRING "${source_line}" ${index} 1 current_char)

      if(in_comment)
        if(NOT current_char STREQUAL "'" AND NOT current_char STREQUAL "\"")
          string(APPEND processed_line "${current_char}")
        endif()
      else()
        if(current_char STREQUAL "'")
          if(in_single_quote)
            math(EXPR next_index "${index} + 1")
            if(next_index LESS line_length)
              string(SUBSTRING "${source_line}" ${next_index} 1 next_char)
              if(next_char STREQUAL "'")
                string(APPEND processed_line "''")
                math(EXPR index "${index} + 2")
                continue()
              endif()
            endif()
            set(in_single_quote FALSE)
          else()
            set(in_single_quote TRUE)
          endif()

          string(APPEND processed_line "'")
        elseif(current_char STREQUAL "!")
          if(NOT in_single_quote)
            set(in_comment TRUE)
          endif()

          string(APPEND processed_line "!")
        else()
          string(APPEND processed_line "${current_char}")
        endif()
      endif()

      math(EXPR index "${index} + 1")
    endwhile()

    string(APPEND output_text "${processed_line}\n")
  endif()
endwhile()

get_filename_component(output_directory "${OUTPUT_FILE}" DIRECTORY)
file(MAKE_DIRECTORY "${output_directory}")
file(WRITE "${OUTPUT_FILE}" "${output_text}")