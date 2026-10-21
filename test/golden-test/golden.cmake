# Runs one golden scenario and compares its output to the reference.
#
# Variables:
#   GOLDEN     path to the golden driver
#   MODE       "values" or "observe"
#   SUITE      suite name
#   OBSERVER   observer name (observe mode)
#   OPTIONS    suite options (observe mode)
#   REFERENCE  reference file
#   WORK_DIR   scratch directory
#   UPDATE     if true, write the reference instead of comparing

file(REMOVE_RECURSE "${WORK_DIR}")
file(MAKE_DIRECTORY "${WORK_DIR}")
set(output "${WORK_DIR}/output.txt")

if(MODE STREQUAL "values")
  execute_process(COMMAND "${GOLDEN}" values "${SUITE}"
                  OUTPUT_FILE "${output}"
                  RESULT_VARIABLE result)
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "golden values ${SUITE} failed: ${result}")
  endif()
elseif(MODE STREQUAL "observe")
  execute_process(COMMAND "${GOLDEN}" observe "${SUITE}" "${OBSERVER}" "${OPTIONS}"
                  WORKING_DIRECTORY "${WORK_DIR}"
                  RESULT_VARIABLE result
                  ERROR_QUIET)
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "golden observe ${SUITE} ${OBSERVER} failed: ${result}")
  endif()
  # One line per logger file with a hash of its normalized content.
  file(GLOB_RECURSE files RELATIVE "${WORK_DIR}/exdata" "${WORK_DIR}/exdata/*")
  list(SORT files)
  set(listing "")
  foreach(file IN LISTS files)
    file(READ "${WORK_DIR}/exdata/${file}" content)
    # The version changes with every commit
    string(REGEX REPLACE "coco_version = '[^']*'" "coco_version = 'VERSION'" content "${content}")
    string(SHA256 digest "${content}")
    string(APPEND listing "${file} ${digest}\n")
  endforeach()
  file(WRITE "${output}" "${listing}")
else()
  message(FATAL_ERROR "Unknown MODE '${MODE}'")
endif()

if(UPDATE)
  get_filename_component(reference_dir "${REFERENCE}" DIRECTORY)
  file(MAKE_DIRECTORY "${reference_dir}")
  file(READ "${output}" content)
  file(WRITE "${REFERENCE}" "${content}")
  message(STATUS "Updated ${REFERENCE}")
  return()
endif()

file(READ "${REFERENCE}" expected)
file(READ "${output}" actual)
if(NOT expected STREQUAL actual)
  string(REPLACE "\n" ";" expected_lines "${expected}")
  string(REPLACE "\n" ";" actual_lines "${actual}")
  set(report "")
  foreach(line IN LISTS actual_lines)
    if(NOT line IN_LIST expected_lines)
      string(APPEND report "  actual:   ${line}\n")
    endif()
  endforeach()
  foreach(line IN LISTS expected_lines)
    if(NOT line IN_LIST actual_lines)
      string(APPEND report "  expected: ${line}\n")
    endif()
  endforeach()
  message(FATAL_ERROR "Output differs from ${REFERENCE}:\n${report}"
                      "Use '${GOLDEN} values ${SUITE} --dump' on both versions to find the difference.")
endif()
