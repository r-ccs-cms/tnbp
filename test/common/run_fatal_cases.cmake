# Abort regressions must fail promptly with our diagnostic, not merely fail to launch.
string(REPLACE "," ";" cases "${CASES}")
foreach(case IN LISTS cases)
    execute_process(COMMAND "${MPIEXEC}" "${NUMPROC_FLAG}" "${RANKS}" "${EXECUTABLE}" fatal "${case}"
        RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error TIMEOUT 20)
    if(NOT "${result}" MATCHES "^[0-9]+$" OR "${result}" STREQUAL "0"
       OR NOT "${output}${error}" MATCHES "tnbp: fatal tensor operation")
        message(FATAL_ERROR "fatal case ${case}: result=${result}\n${output}\n${error}")
    endif()
endforeach()
message(STATUS "Fatal cases ${CASES} terminate promptly with diagnostics (${RANKS} ranks)")
