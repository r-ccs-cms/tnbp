# Validate application wiring, without stochastic updates obscuring the input:
# two singleton lines are skipped by two-site sampling in both layout modes.
file(MAKE_DIRECTORY "${DATA}")
file(WRITE "${DATA}/edges" "7 20\n")
file(WRITE "${DATA}/lines" "7\n20\n")
math(EXPR owner "${RANKS}-1")
file(WRITE "${DATA}/owners" "7 ${owner}\n20 ${owner}\n")
set(base "${MPIEXEC}" "${NUMPROC_FLAG}" "${RANKS}" "${EXECUTABLE}"
    --inputedge "${DATA}/edges" --linefile "${DATA}/lines" --ownerfile "${DATA}/owners"
    --initial_statefile "${DATA}/states" --num_batches 2 --update_sites 2 --p_bond_dim 16)
file(WRITE "${DATA}/states" "# site batch0 batch1\n20 5 6 # comment\n\n7 3 12\n")
foreach(mode IN ITEMS single alternating)
    set(extra)
    if(mode STREQUAL "alternating")
        set(extra --second_linefile "${DATA}/lines" --first_round_trips 2 --second_round_trips 3)
    endif()
    execute_process(COMMAND ${base} --line_mode ${mode} ${extra} --num_samples 2
        RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error TIMEOUT 30)
    if(NOT "${result}" STREQUAL "0" OR NOT output MATCHES "sample 0.0: 35"
       OR NOT output MATCHES "sample 0.1: 126" OR NOT output MATCHES "sample 1.1: 126")
        message(FATAL_ERROR "initial states ${mode}: ${result}\n${output}\n${error}")
    endif()
endforeach()
# Validation still runs with zero samples. Root parsing and owner dimension errors
# must terminate the whole job promptly, including when rank zero has no tensors.
set(contents "7 3 12\n7 5 6\n" "7 3 12\n" "7 3 12\n21 5 6\n"
    "7 3\n20 5 6\n" "7 3 12 1\n20 5 6\n" "7 -1 12\n20 5 6\n"
    "7 3junk 12\n20 5 6\n" "7 999999999999999999 12\n20 5 6\n" "7 3 16\n20 5 6\n")
foreach(content IN LISTS contents)
    file(WRITE "${DATA}/states" "${content}")
    execute_process(COMMAND ${base} --num_samples 0 RESULT_VARIABLE result
        OUTPUT_VARIABLE output ERROR_VARIABLE error TIMEOUT 20)
    if(NOT "${result}" MATCHES "^[0-9]+$" OR "${result}" STREQUAL "0"
       OR NOT "${error}" MATCHES "initial state file:")
        message(FATAL_ERROR "invalid initial states: ${result}\n${output}\n${error}")
    endif()
endforeach()
execute_process(COMMAND ${base} --initial_statefile "${DATA}/absent" --num_samples 0
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error TIMEOUT 20)
if(NOT "${result}" MATCHES "^[0-9]+$" OR "${result}" STREQUAL "0"
   OR NOT error MATCHES "cannot open initial state file")
    message(FATAL_ERROR "missing initial states: ${result}\n${output}\n${error}")
endif()
message(STATUS "initial state CLI PASS (${RANKS} ranks): both modes, two batches, 16 states, input errors")
