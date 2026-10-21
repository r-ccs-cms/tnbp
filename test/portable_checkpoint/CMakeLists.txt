add_executable(estimator_portable ../../apps/estimator_restart/main.cc)
target_link_libraries(estimator_portable PRIVATE tnbp_test_backend)
if(TNBP_TEST_BACKEND STREQUAL "cuda")
  target_compile_definitions(estimator_portable PRIVATE USE_TCAPI_CUDA)
endif()
add_executable(portable_checkpoint main.cc)
target_compile_definitions(portable_checkpoint PRIVATE USE_COMPLEX)
target_link_libraries(portable_checkpoint PRIVATE tnbp_test_backend)
if(TNBP_TEST_BACKEND STREQUAL "cuda")
  target_compile_definitions(portable_checkpoint PRIVATE USE_TCAPI_CUDA)
endif()
if(BUILD_TESTING AND TNBP_TEST_BACKEND STREQUAL "cpu")
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
  add_test(NAME portable_checkpoint_acceptance COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/integration.py ${MPIEXEC_EXECUTABLE}
    $<TARGET_FILE:portable_checkpoint> $<TARGET_FILE:estimator_portable>
    $<TARGET_FILE:mcmc_sampling_complex> $<TARGET_FILE:mcmc_sampling>
    ${CMAKE_CURRENT_BINARY_DIR}/portable-evidence)
  set_tests_properties(portable_checkpoint_acceptance PROPERTIES TIMEOUT 600 PROCESSORS 4)
endif()
