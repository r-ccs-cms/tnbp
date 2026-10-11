cmake_minimum_required(VERSION 3.18)
# Retain the previous configure entry point for existing build directories.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(tnbp_backend_integration LANGUAGES CXX)
    include(CTest)
    include(../cmake/Backend.cmake)
endif()
add_executable(backend_integration main.cc)
target_link_libraries(backend_integration PRIVATE tnbp_test_backend)
# CPU tests use CMake's MPI launcher. CUDA tests must be launched explicitly
# inside a scheduler allocation with one visible GPU per process.
if(BUILD_TESTING AND TNBP_TEST_BACKEND STREQUAL "cpu")
    foreach(ranks IN ITEMS 1 2)
        add_test(NAME integration_${ranks}_rank COMMAND ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:backend_integration> ${MPIEXEC_POSTFLAGS})
        set_tests_properties(integration_${ranks}_rank PROPERTIES TIMEOUT 120 PROCESSORS ${ranks}
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
    endforeach()
endif()

# Compile each public TNBP header as the first and only include in its own
# translation unit. Also check all headers in reverse order. Object files are
# deliberately not linked: this is an include-dependency check, not an ODR test.
file(GLOB_RECURSE TNBP_PUBLIC_HEADERS CONFIGURE_DEPENDS
    RELATIVE "${TNBP_ROOT}/include" "${TNBP_ROOT}/include/tnbp/*.h"
    "${TNBP_ROOT}/include/qasm/*.h" "${TNBP_ROOT}/include/pauli/*.h")
set(header_check_sources)
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/header_checks")
foreach(header IN LISTS TNBP_PUBLIC_HEADERS)
    string(REPLACE "/" "_" unit "${header}")
    set(source "${CMAKE_CURRENT_BINARY_DIR}/header_checks/${unit}.cc")
    file(WRITE "${source}" "#include \"${header}\"\n")
    list(APPEND header_check_sources "${source}")
endforeach()
list(REVERSE TNBP_PUBLIC_HEADERS)
set(reverse_source "${CMAKE_CURRENT_BINARY_DIR}/header_checks/reverse_order.cc")
file(WRITE "${reverse_source}" "")
foreach(header IN LISTS TNBP_PUBLIC_HEADERS)
    file(APPEND "${reverse_source}" "#include \"${header}\"\n")
endforeach()
list(APPEND header_check_sources "${reverse_source}")
add_library(header_checks OBJECT ${header_check_sources})
target_link_libraries(header_checks PRIVATE tnbp_test_backend)

# Link three independently compiled users of the public headers. Do not allow
# optimization/unity builds to hide ordinary external definitions.
add_executable(header_link_check link_main.cc link_a.cc link_b.cc)
target_link_libraries(header_link_check PRIVATE tnbp_test_backend)
set_target_properties(header_link_check PROPERTIES UNITY_BUILD OFF)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(header_link_check PRIVATE -O0 -fno-inline)
endif()
if(BUILD_TESTING)
    add_test(NAME header_link_check COMMAND header_link_check)
endif()


if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    add_subdirectory(../graph_tensor_contraction suites/graph_tensor_contraction)
    add_subdirectory(../boundary_mps suites/boundary_mps)
    add_subdirectory(../mcmc_sampling suites/mcmc_sampling)
endif()
