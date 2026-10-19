add_executable(bp_trace_phase main.cc)
target_link_libraries(bp_trace_phase PRIVATE tnbp_test_backend)
if(BUILD_TESTING AND TNBP_TEST_BACKEND STREQUAL "cpu")
  foreach(ranks IN ITEMS 1 2 3 4)
    add_test(NAME bp_trace_phase_${ranks}_rank COMMAND ${MPIEXEC_EXECUTABLE}
      ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:bp_trace_phase> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(bp_trace_phase_${ranks}_rank PROPERTIES TIMEOUT 120 PROCESSORS ${ranks}
      ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
  endforeach()
endif()

add_executable(bp_message_difference message_difference.cc)
target_link_libraries(bp_message_difference PRIVATE tnbp_test_backend)
if(BUILD_TESTING AND TNBP_TEST_BACKEND STREQUAL "cpu")
  foreach(ranks IN ITEMS 1 2 4)
    add_test(NAME bp_message_difference_${ranks}_rank COMMAND ${MPIEXEC_EXECUTABLE}
      ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:bp_message_difference> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(bp_message_difference_${ranks}_rank PROPERTIES TIMEOUT 120 PROCESSORS ${ranks}
      ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
  endforeach()
endif()
