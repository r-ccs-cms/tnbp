cmake_minimum_required(VERSION 3.18)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(tnbp_pauli_evolution LANGUAGES CXX)
  include(CTest)
  include(../cmake/Backend.cmake)
endif()
add_executable(pauli_evolution main.cc)
target_link_libraries(pauli_evolution PRIVATE tnbp_test_backend)
if(BUILD_TESTING AND TNBP_TEST_BACKEND STREQUAL "cpu")
  foreach(ranks IN ITEMS 1 2)
    add_test(NAME pauli_evolution_${ranks}_rank COMMAND ${MPIEXEC_EXECUTABLE}
      ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:pauli_evolution> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(pauli_evolution_${ranks}_rank PROPERTIES TIMEOUT 120 PROCESSORS ${ranks}
      ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
  endforeach()
endif()
add_subdirectory(../../apps/imaginary_time_evolution imaginary_time_cli)
if(BUILD_TESTING AND TNBP_TEST_BACKEND STREQUAL "cpu")
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
  add_test(NAME imaginary_time_cli COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/cli.py ${MPIEXEC_EXECUTABLE}
    $<TARGET_FILE:imaginary_time_evolution> ${CMAKE_CURRENT_BINARY_DIR}/cli-evidence)
  set_tests_properties(imaginary_time_cli PROPERTIES TIMEOUT 240 PROCESSORS 2)
endif()
add_executable(truncation_batch truncation_batch.cc)
target_link_libraries(truncation_batch PRIVATE tnbp_test_backend)
if(BUILD_TESTING AND TNBP_TEST_BACKEND STREQUAL "cpu")
  foreach(ranks IN ITEMS 1 2 3 4)
    add_test(NAME truncation_batch_${ranks}_rank COMMAND ${MPIEXEC_EXECUTABLE}
      ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:truncation_batch> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(truncation_batch_${ranks}_rank PROPERTIES TIMEOUT 120 PROCESSORS ${ranks}
      ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
  endforeach()
endif()
add_executable(graph_input graph_input.cc)
target_link_libraries(graph_input PRIVATE tnbp_test_backend)
if(BUILD_TESTING AND TNBP_TEST_BACKEND STREQUAL "cpu")
  add_test(NAME graph_input COMMAND graph_input ${CMAKE_CURRENT_BINARY_DIR})
endif()
