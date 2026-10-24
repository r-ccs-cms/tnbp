add_executable(checkpoint_file_paths file_paths.cc)
target_include_directories(checkpoint_file_paths PRIVATE "${TNBP_ROOT}/include")
target_compile_features(checkpoint_file_paths PRIVATE cxx_std_17)
if(BUILD_TESTING)
  add_test(NAME checkpoint_file_paths COMMAND checkpoint_file_paths)
endif()
add_executable(estimator_checkpoint_fixture estimator_checkpoint.cc)
target_compile_definitions(estimator_checkpoint_fixture PRIVATE USE_COMPLEX)
target_link_libraries(estimator_checkpoint_fixture PRIVATE tnbp_test_backend)
if(TNBP_TEST_BACKEND STREQUAL "cuda")
  target_compile_definitions(estimator_checkpoint_fixture PRIVATE USE_TCAPI_CUDA)
endif()
add_executable(estimator_cli ../../apps/estimator/main.cc)
target_link_libraries(estimator_cli PRIVATE tnbp_test_backend)
if(TNBP_TEST_BACKEND STREQUAL "cuda")
  target_compile_definitions(estimator_cli PRIVATE USE_TCAPI_CUDA)
endif()
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
  add_test(NAME checkpoint_review_regression COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/review_regression.py ${MPIEXEC_EXECUTABLE}
    $<TARGET_FILE:estimator_portable> $<TARGET_FILE:estimator_checkpoint_fixture>
    $<TARGET_FILE:portable_checkpoint> ${CMAKE_CURRENT_BINARY_DIR}/review-evidence)
  set_tests_properties(checkpoint_review_regression PROPERTIES TIMEOUT 300 PROCESSORS 2)
  add_test(NAME portable_checkpoint_acceptance COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/integration.py ${MPIEXEC_EXECUTABLE}
    $<TARGET_FILE:portable_checkpoint> $<TARGET_FILE:estimator_portable>
    $<TARGET_FILE:mcmc_sampling_complex> $<TARGET_FILE:mcmc_sampling>
    ${CMAKE_CURRENT_BINARY_DIR}/portable-evidence)
  set_tests_properties(portable_checkpoint_acceptance PROPERTIES TIMEOUT 600 PROCESSORS 4)
  add_test(NAME app_site_ownership COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/ownership.py ${MPIEXEC_EXECUTABLE}
    $<TARGET_FILE:estimator_cli> $<TARGET_FILE:estimator_portable>
    $<TARGET_FILE:mcmc_sampling_complex> $<TARGET_FILE:mcmc_sampling>
    ${CMAKE_CURRENT_BINARY_DIR}/ownership-evidence)
  set_tests_properties(app_site_ownership PROPERTIES TIMEOUT 300 PROCESSORS 4)
endif()
