# Full-SVD relocation of long boundary bonds, with optional CUDA compilation.
add_executable(bmps_shift shift_regression.cc)
target_link_libraries(bmps_shift PRIVATE tnbp_test_backend)
if(BUILD_TESTING AND TNBP_TEST_BACKEND STREQUAL "cpu")
    foreach(ranks IN ITEMS 1 2 3)
        add_test(NAME bmps_shift_${ranks}_rank COMMAND ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:bmps_shift> ${MPIEXEC_POSTFLAGS})
        set_tests_properties(bmps_shift_${ranks}_rank PROPERTIES TIMEOUT 90 PROCESSORS ${ranks}
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
    endforeach()
endif()

if(TNBP_TEST_BACKEND STREQUAL "cpu")
    add_executable(bmps_upper_independence upper_independence.cc)
    target_link_libraries(bmps_upper_independence PRIVATE tnbp_test_backend)
    if(BUILD_TESTING)
        foreach(ranks IN ITEMS 1 2)
            add_test(NAME bmps_upper_independence_${ranks}_rank COMMAND ${MPIEXEC_EXECUTABLE}
                ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:bmps_upper_independence> ${MPIEXEC_POSTFLAGS})
            set_tests_properties(bmps_upper_independence_${ranks}_rank PROPERTIES TIMEOUT 60 PROCESSORS ${ranks}
                ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
        endforeach()
    endif()
endif()
