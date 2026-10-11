# Generic graph contraction and bond relocation, with optional CUDA compilation.
add_executable(graph_contraction_regression regression.cc)
target_link_libraries(graph_contraction_regression PRIVATE tnbp_test_backend)
if(BUILD_TESTING AND TNBP_TEST_BACKEND STREQUAL "cpu")
    foreach(ranks IN ITEMS 1 2 3)
        add_test(NAME graph_contraction_regression_${ranks}_rank COMMAND ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:graph_contraction_regression> ${MPIEXEC_POSTFLAGS})
        set_tests_properties(graph_contraction_regression_${ranks}_rank PROPERTIES TIMEOUT 90 PROCESSORS ${ranks}
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
    endforeach()
endif()
